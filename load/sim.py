# -*- coding: utf-8 -*-
"""
Created on Tue Jan 13 23:27:16 2015

@author: hoseung
"""

# coding: utf-8
import numpy as np

class Simbase:
    """ base    
    """
    def __init__(self):
        
        pass
    
    def add_amr(self):
        from load import amr
        self.amr = amr.Amr(self.info)
        print("An AMR instance is created\n")

    def _hilbert_cpulist(self, info, ranges):
        '''
        After determining cpu numbers, read the cpu files and cut off data
        that are outside ranges.
        -> cpu files contain data points within a given ragnes
        BUT! they also contain other data points outside the ragnes.
        '''
        if not(hasattr(self, 'amr')):
            print("No AMR instance,")
            print("Loading one...")
            self.add_amr()
        nlevelmax = self.amr.header.nlevelmax
        nboundary = self.amr.header.nboundary
        ndim = self.info.ndim
        ncpu = self.info.ncpu_tot  #ncpu_tot
        # header of part output is needed (and is prepared.)

        # header of amr output is needed.
        two2ndim = 2**ndim
        nx = self.amr.header.ng[0]
        ny = self.amr.header.ng[1]
        nz = self.amr.header.ng[2]

        xbound = [nx/2., ny/2., nz/2.]
        ngridfile = np.zeros(ncpu + nboundary, nlevelmax)
        ngridlevel = np.zeros(ncpu + nlevelmax)

        if (nboundary > 0):
            ngridbound = np.zeros(nboundary, nlevelmax)
        # info file is needed. (and is ready)
        # cpu_list & hilbert_key : sim.info.

        try:
            lmax
        except:
            lmax = nlevelmax
#        print(' >>> working resolution (lmax) =', lmax)

        xxmin = ranges[0][0]
        xxmax = ranges[0][1]
        yymin = ranges[1][0]
        yymax = ranges[1][1]
        zzmin = ranges[2][0]
        zzmax = ranges[2][1]

        dmax = max([xxmax-xxmin, yymax-yymin, zzmax-zzmin])
        for ilevel in range(1, lmax):
            dx = 0.5**ilevel
            if (dx <= dmax):
                lmin = ilevel
                break

        bit_length = lmin - 1
        maxdom = 2**bit_length
        imin = 0
        imax = 0
        jmin = 0
        jmax = 0
        kmin = 0
        kmax = 0

        if (bit_length >= 0):
            imin = int(xxmin * maxdom)
            imax = imin + 1
            jmin = int(yymin * maxdom)
            jmax = jmin + 1
            kmin = int(zzmin * maxdom)
            kmax = kmin + 1

        dkey = (2**(nlevelmax+1)/maxdom)**(ndim)
        ndom = 1
        if (bit_length > 0):
            ndom = 8
        idom = np.zeros(9)
        jdom = np.zeros(9)
        kdom = np.zeros(9)
        idom[1] = imin; idom[2] = imax; idom[3] = imin; idom[4] = imax
        idom[5] = imin; idom[6] = imax; idom[7] = imin; idom[8] = imax
        jdom[1] = jmin; jdom[2] = jmin; jdom[3] = jmax; jdom[4] = jmax
        jdom[5] = jmin; jdom[6] = jmin; jdom[7] = jmax; jdom[8] = jmax
        kdom[1] = kmin; kdom[2] = kmin; kdom[3] = kmin; kdom[4] = kmin
        kdom[5] = kmax; kdom[6] = kmax; kdom[7] = kmax; kdom[8] = kmax

        bounding_min = np.zeros(9)
        bounding_max = np.zeros(9)

        for i in range(1, ndom + 1):
            if bit_length > 0:
                order_min = self._hilbert3d([idom[i]], [jdom[i]], [kdom[i]],
                                            bit_length, 1)
                # order_min, array or single variable??
                # Will it be vectorized??
            else:
                order_min = 0

            order_min = np.asarray(order_min)
            bounding_min[i] = [order_min][0] * dkey
            bounding_max[i] = ([order_min][0] + 1) * dkey
            # [(x + 1) * dkey for x in order_min]

        cpu_min = np.zeros(9, dtype=np.int)
        cpu_max = np.zeros(9, dtype=np.int)
        bound_key = info.hilbertkey[0]

        cpu_list = np.zeros(ncpu+1, dtype=np.int)

        for impi in range(1, ncpu + 1):
            for i in range(1, ndom + 1):
                if ((bound_key[impi - 1] <= bounding_min[i]) and
                   (bound_key[impi] > bounding_min[i])):
                    cpu_min[i] = impi
                if ((bound_key[impi - 1] < bounding_max[i]) and
                   (bound_key[impi] >= bounding_max[i])):
                    cpu_max[i] = impi

        ncpu_read = 0
        cpu_read = np.zeros(ncpu + 1, dtype=np.int)
        for i in range(1, ndom + 1):
            for j in range(cpu_min[i], cpu_max[i] + 1):
                # np.arange(10,10) = [], np.arange(10,11) = [10]
                if cpu_read[j] == 0:
                    ncpu_read += 1
                    cpu_list[ncpu_read] = j
                    cpu_read[j] = 1

        for ilevel in range(1, lmax + 1):
            nx_full = 2**ilevel
            ny_full = 2**ilevel
            nz_full = 2**ilevel
            imin = xxmin * (nx_full) + 1
            imax = xxmax * (nx_full) + 1
            jmin = yymin * (ny_full) + 1
            jmax = yymax * (ny_full) + 1
            kmin = zzmin * (nz_full) + 1
            kmax = zzmax * (nz_full) + 1

# Sort cpu_list in descending npart order for memory efficiency
        cpu_list = cpu_list[cpu_list > 0]  # crop empty part

        return np.sort(cpu_list)

    def _hilbert3d(self, x, y, z, bit_length, npoint):
        '''
        Calculate hilbert doamin decomposition
        '''
        state_diagram = np.zeros((8, 2, 12), dtype=np.int)
        state_diagram[:, 0, 0] = [1, 2, 3, 2, 4, 5, 3, 5]
        state_diagram[:, 1, 0] = [0, 1, 3, 2, 7, 6, 4, 5]
        state_diagram[:, 0, 1] = [2, 6, 0, 7, 8, 8, 0, 7]
        state_diagram[:, 1, 1] = [0, 7, 1, 6, 3, 4, 2, 5]
        state_diagram[:, 0, 2] = [0, 9, 10, 9, 1, 1, 11, 11]
        state_diagram[:, 1, 2] = [0, 3, 7, 4, 1, 2, 6, 5]
        state_diagram[:, 0, 3] = [6, 0, 6, 11, 9, 0, 9, 8]
        state_diagram[:, 1, 3] = [2, 3, 1, 0, 5, 4, 6, 7]
        state_diagram[:, 0, 4] = [11, 11, 0, 7, 5, 9, 0, 7]
        state_diagram[:, 1, 4] = [4, 3, 5, 2, 7, 0, 6, 1]
        state_diagram[:, 0, 5] = [4, 4, 8, 8, 0, 6, 10, 6]
        state_diagram[:, 1, 5] = [6, 5, 1, 2, 7, 4, 0, 3]
        state_diagram[:, 0, 6] = [5, 7, 5, 3, 1, 1, 11, 11]
        state_diagram[:, 1, 6] = [4, 7, 3, 0, 5, 6, 2, 1]
        state_diagram[:, 0, 7] = [6, 1, 6, 10, 9, 4, 9, 10]
        state_diagram[:, 1, 7] = [6, 7, 5, 4, 1, 0, 2, 3]
        state_diagram[:, 0, 8] = [10, 3, 1, 1, 10, 3, 5, 9]
        state_diagram[:, 1, 8] = [2, 5, 3, 4, 1, 6, 0, 7]
        state_diagram[:, 0, 9] = [4, 4, 8, 8, 2, 7, 2, 3]
        state_diagram[:, 1, 9] = [2, 1, 5, 6, 3, 0, 4, 7]
        state_diagram[:, 0, 10] = [7, 2, 11, 2, 7, 5, 8, 5]
        state_diagram[:, 1, 10] = [4, 5, 7, 6, 3, 2, 0, 1]
        state_diagram[:, 0, 11] = [10, 3, 2, 6, 10, 3, 4, 4]
        state_diagram[:, 1, 11] = [6, 1, 7, 0, 5, 2, 4, 3]

        i_bit_mask = np.zeros(3 * bit_length)
        ind = np.arange(bit_length)
#       order      = dblarr(npoint)

        for ip in range(npoint):  # check if loop range is valid.
            # convert to binary
            x_bit_mask = self._btest(x[ip], bit_length - 1, True)
            y_bit_mask = self._btest(y[ip], bit_length - 1, True)
            z_bit_mask = self._btest(z[ip], bit_length - 1, True)

            # interleave bits
            i_bit_mask[3 * ind + 2] = x_bit_mask
            i_bit_mask[3 * ind + 1] = y_bit_mask
            i_bit_mask[3 * ind + 0] = z_bit_mask

            # build Hilbert ordering using state diagram
            cstate = 0
            # from bit_length -1 to 0 in descending order.
            for i in range(bit_length - 1, -1, -1):
                b2 = 0
                if (i_bit_mask[3 * i + 2]):
                    b2 = 1
                b1 = 0
                if (i_bit_mask[3 * i + 1]):
                    b1 = 1
                b0 = 0
                if (i_bit_mask[3 * i + 0]):
                    b0 = 1
                sdigit = b2 * 4 + b1 * 2 + b0

                nstate = state_diagram[sdigit, 0, cstate]
                hdigit = state_diagram[sdigit, 1, cstate]
                i_bit_mask[3 * i + 2] = self._btest(hdigit, 2, p_all=False)
                i_bit_mask[3 * i + 1] = self._btest(hdigit, 1, p_all=False)
                i_bit_mask[3 * i + 0] = self._btest(hdigit, 0, p_all=False)
                cstate = nstate

            # save Hilbert key as double precision real
            order = [0.0] * npoint
            for i in range(3 * bit_length):
                b0 = 0
                if (i_bit_mask[i]):
                    b0 = 1
                order[ip] += b0 * (2**i)
        return order

    def _btest(self, tmp, bit, p_all=False):
        nbit = bit
        if (not p_all and tmp != 0):
            tmp2 = int(np.log2(tmp))+1
            if (tmp2 > nbit):
                nbit = tmp2

        res = [0]*(nbit+1)

        for j in np.arange(nbit + 1, 0, -1) - 1:
            res[j] = np.int(tmp / 2**j)
            tmp -= res[j] * 2**j

        if (p_all):
            return(res)
        else:
            return(res[bit])

#%%
# Sim class has (logically linked) classes: Info, Part, Amr
class Sim(Simbase):
    """
    Defines the 'host class' of part, amr, hydro, info.

    Global information is stored in this class: 
        ndim, ncpu, base, type of simulation (DMO, zoom, and so on)
    Later it will also include .nml information.
    (Romain's git version generates such output in text files)

    Currently this class deals with single snapshot.
    But I hope to expand it for multiple snapshots.

    .. note::
        This is a sample note
    """
    def __init__(self, nout=None, base='./', data_dir='snapshots',
                 ranges=[[0.0,1.0],[0.0,1.0],[0.0,1.0]], dmo=False
                 , setup=False, region=None):

        self.nout = nout            
        self.set_base(base)
        # info appreciates nout and base (not medatary, though)
        self.dmo = dmo
        # DMO runs are slightly different!
        self.add_info()
        self.set_data_dir(data_dir)
        # set_data_dir and set_range needs info instance be exist.
#        self.set_ranges(ranges)
        if region is not None:
            ranges = [[region['xc'] - region['radius'], region['xc'] + region['radius']],
                      [region['yc'] - region['radius'], region['yc'] + region['radius']],
                      [region['zc'] - region['radius'], region['zc'] + region['radius']]]
        if setup:
            self.setup(nout, base, data_dir, ranges, dmo)

    def __call__(self, *args):
        # Function emulation      
        return self.__init__(*args)
        self.setup(self, *args)

    def all_set(self):
        return (self.nout is not None) & (self.base is not None)

    def setup(self, nout=None, base='./', data_dir='snapshots',
                 ranges=[[0.0,1.0],[0.0,1.0],[0.0,1.0]], dmo=False):
        self.nout = nout
        self.set_base(base)        
        if self.nout is None:
            print("Note that 'nout' is not set. \n use sim.Sim.set_nout(nout)")        
        self.add_info()
        # set_data_dir and set_range needs info instance be exist.
        self.set_ranges(ranges)

        if self.all_set():
            self.add_amr()
            self.set_cpus(self._hilbert_cpulist(self.info, self.ranges))
        print(' Simulation set up.')

    def set_base(self, base):
        """
            Sets Working directory.
        """
        from os import path
        self.base = path.abspath(base)
        #self.show_base()

    def set_nout(self, nout):
        """
            Sets output number.
            a list of nouts will be supported in the future, soon, I hope.
        """
        self.nout = nout

    def set_data_dir(self, data_dir):
        """
        By default, row data is in simulation_base/snapshots/
        """
        from os import path
        self.data_dir =  path.join(self.base, '', data_dir, '')
#        self.info.set_data_dir(self.data_dir)

    def show_base(self):
        print("setting the base(working) directory to :", self.base)

    def add_hydro(self, load=False, lmax=19):
        from load import hydro
        self.hydro = hydro.Hydro(self.info, self.amr)
        print("An Hydro instance is created\n")
        if load :
            self.hydro.amr2cell(lmax=lmax)
        else:
            print("Use hydro.amr2cell() to load hydro variables")
            
    def add_info(self, load=False):
        from load import info
        self.info = info.Info(self.nout, self.base, load=load)
#        self.info.setup()

    def add_part(self, ptypes=[], load=False, fortran=True, dmo=False, **kwargs):
        """
        Add particle instance to the simulation instance. 
        Requires types of particles and particle data.
        load = True  to load actual data on creating the instance

        parameters
        ----------
        ptypes : list of particle type and information. 
                ["dm id pos"] or ["dm id pos", "star mass vel"]
                
        """
        if dmo:
            self.dmo = True
        from load import part
        print("Types of particles you want to load are: ", ptypes)
        self.part = part.Part(self.info, ptypes=ptypes, data_dir=self.data_dir, dmo=self.dmo, **kwargs)
        print("A particle instance is created\n")

        if load :
            self.part.load(fortran=fortran)
        else:
            print("Use part.load() to load particle")

    def get_cpus(self):
        return self.cpus

    def set_cpus(self, cpus):
        self.cpus = cpus
        try:
            print("Updating info.cpus")
            self.info._set_cpus(self.get_cpus())
        except AttributeError:
            print("No info._set_cpus attribute??")
        #self.show_cpus()

    def show_cpus(self):
        print(" ncpus : %s \n" % self.get_cpus())

    def set_ranges(self, ranges=[[0, 1], [0, 1], [0, 1]]):
        nr = np.asarray(ranges)
        if not(nr.shape[0] == 3 and nr.shape[1] == 2):
            # Because actual operation on the given input(ranges)
            # does not take place soon, it's not a good place to use
            # try & except clause. There is nothing to try yet.
            print(' Error!')
            print('Shape of ranges is wrong:', nr.shape)
            print('example : [[0.1,0.3],[0.2,0.4],[0.6,0.8]] \n')
        else:
            self.ranges = ranges
            try:
                self.info._set_ranges(self.ranges)
            except AttributeError:
                print("There is no info._set_ranges attribute")
            self.show_ranges()
            self.set_cpus(self._hilbert_cpulist(self.info, self.ranges))

    def show_ranges(self):
        print("Ranges = %s\n" % self.ranges)

    def search_zoomin_region(self, *args, **kwargs):  # If part is not loaded yet, load particles
        """
        Not only part, but also hydro or amr can be used to find the zoomin region!

        Determine priority, amr or part?
        """
        if hasattr(self, 'part'):
            print("Have part")
            self.set_zregion(self.part.search_zoomin( *args, **kwargs))

        if hasattr(self, 'amr'):
            print("Have amr")
            self.set_zregion(self.amr.search_zoomin( *args, **kwargs))

    def set_zregion(self, zregion):
        self.zregion = zregion
        self.info.zregion = zregion



