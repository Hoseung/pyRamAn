# -*- coding: utf-8 -*-
"""
Created on Tue Jan 13 23:27:16 2015

@author: hoseung
"""

# coding: utf-8
import numpy as np

# Sim class has (logically linked) classes: Info, Part, Amr
#
#


class Sim(object):

    def __init__(self, nout, base='./', ranges=[[0.0,1.0],[0.0,1.0],[0.0,1.0]]):
        self.nout = nout
        self.base = base
        self.info = Info(self.nout, self.base)
        self.ranges = ranges
        self.cpus = self._hilbert_cpulist(self.info, self.ranges)
        #print('cpu list:', self.cpus)

    def add_part(self, ptypes):
        print("Types of particles you want to load are: ", ptypes)
        self.part = Part(self.info, ptypes)
        print("A particle instance is created")
        print("Use part.load() to load these particle")

    def add_amr(self):
        self.amr = Amr(self.info)
        print("An AMR instance is created")

    def add_hydro(self):
        self.hydro = Hydro(self.info)
        print("An Hydro instance is created")

    def set_ranges(self, new_ranges=[[0,1],[0,1],[0,1]]):
        nr = np.asarray(new_ranges)
        if not(nr.shape[0] == 3 and nr.shape[1] == 2):
            print('Shape of ranges is wrong:', nr.shape)
            print('example : [[0.1,0.3],[0.2,0.4],[0.6,0.8]]')
        else:
            self.ranges = new_ranges
        print(self.ranges)
        self.cpus = self._hilbert_cpulist(self.info, new_ranges)

    def zoomin_region(self):        # If part is not loaded yet, load particles
        self.part.get_zoomin()
        
        #ind = np.where(mass == min_mass)

    def _hilbert_cpulist(self, info, ranges):
        '''
        After determining cpu numbers, read the cpu files and cut off data
        that are outside ranges.
        -> cpu files contain data points within a given ragnes
        BUT! they also contain other data points outside the ragnes.
        '''
        if not(hasattr(self, 'amr')):
            self.add_amr()
        nlevelmax = self.amr.header.nlevelmax
        nboundary = self.amr.header.nboundary
        ndim = self.info.ndim
        ncpu = self.info.ncpu
        # header of part output is needed (and is prepared.)

        # header of amr output is needed.
        two2ndim = 2**ndim
        nx = self.amr.header.ng[0]
        ny = self.amr.header.ng[1]
        nz = self.amr.header.ng[2]

        xbound=[nx/2., ny/2., nz/2.]
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
        print(' >>> working resolution (lmax) =', lmax)

        xxmin = ranges[0][0]
        xxmax = ranges[0][1]
        yymin = ranges[1][0]
        yymax = ranges[1][1]
        zzmin = ranges[2][0]
        zzmax = ranges[2][1]

        dmax=max([xxmax-xxmin, yymax-yymin, zzmax-zzmin])
        for ilevel in range(1, lmax):
            dx = 0.5**ilevel
            if (dx <= dmax):
                lmin = ilevel
                break
        

        bit_length = lmin-1
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
            jmin = int(yymin *maxdom)
            jmax = jmin + 1
            kmin = int(zzmin * maxdom)
            kmax = kmin +1

        dkey = (2**(nlevelmax+1)/maxdom)**(ndim)
        ndom = 1
        if (bit_length > 0): 
            ndom=8
        idom = np.zeros(9)
        jdom = np.zeros(9)
        kdom = np.zeros(9)
        idom[1] = imin ; idom[2] = imax ; idom[3] = imin ; idom[4] = imax
        idom[5] = imin ; idom[6] = imax ; idom[7] = imin ; idom[8] = imax
        jdom[1] = jmin ; jdom[2] = jmin ; jdom[3] = jmax ; jdom[4] = jmax
        jdom[5] = jmin ; jdom[6] = jmin ; jdom[7] = jmax ; jdom[8] = jmax
        kdom[1] = kmin ; kdom[2] = kmin ; kdom[3] = kmin ; kdom[4] = kmin
        kdom[5] = kmax ; kdom[6] = kmax ; kdom[7] = kmax ; kdom[8] = kmax

        bounding_min = np.zeros(9)
        bounding_max = np.zeros(9)
        #print('bit_length:',bit_length)

        for i in range(1,ndom +1):
            if bit_length > 0:
                order_min = self._hilbert3d( [idom[i]], [jdom[i]], [kdom[i]], bit_length, 1)
                # order_min, array or single variable?? 
                # Someday will it be vectorized??               
            else:
                order_min = 0
            #help(dkey) # float
            #help(order_min) # list
            #help(bounding_min) # ndarray
            #bounding_min[i] = [x * dkey for x in order_min]
            order_min = np.asarray(order_min)
            #print(order_min, dkey)
            bounding_min[i] = [order_min][0] * dkey
            bounding_max[i] = ([order_min][0] + 1) * dkey
            #[(x + 1) * dkey for x in order_min] 
        

        cpu_min = np.zeros(9, dtype=np.int)
        cpu_max = np.zeros(9, dtype=np.int)
        bound_key = info.hilbertkey[0]

        cpu_list = np.zeros(ncpu+1, dtype=np.int) # integer
        
        for impi in range(1, ncpu + 1):
            for i in range(1, ndom + 1):
                if (bound_key[impi - 1] <= bounding_min[i]) and (bound_key[impi] > bounding_min[i]):
                    cpu_min[i]=impi
                if (bound_key[impi - 1] < bounding_max[i]) and (bound_key[impi] >= bounding_max[i]):
                    cpu_max[i]=impi
                    
        ncpu_read = 0
        cpu_read  = np.zeros(ncpu + 1, dtype=np.int)
        for i in range(1, ndom + 1):
            for j in range(cpu_min[i], cpu_max[i] + 1): # np.arange(10,10) = [], np.arange(10,11) = [10]
                if cpu_read[j] == 0:
                    ncpu_read += 1
                    cpu_list[ncpu_read] = j
                    cpu_read[j] = 1
        
        for ilevel in range(1, lmax + 1):
            nx_full = 2**ilevel
            ny_full = 2**ilevel
            nz_full = 2**ilevel
            imin    = xxmin * (nx_full) + 1
            imax    = xxmax * (nx_full) + 1
            jmin    = yymin * (ny_full) + 1
            jmax    = yymax * (ny_full) + 1
            kmin    = zzmin * (nz_full) + 1
            kmax    = zzmax * (nz_full) + 1
            
# Sort cpu_list in descending npart order for memory efficiency
        cpu_list = cpu_list[cpu_list > 0] # crop empty part
        
        return np.sort(cpu_list)#[::-1]

    def _hilbert3d(self, x, y, z, bit_length, npoint):
        '''
        Calculate hilbert doamin decomposition 
        '''

        state_diagram=np.zeros((8,2,12),dtype=np.int)
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
        state_diagram[:, 0, 11] = [10,3, 2, 6, 10, 3, 4, 4]
        state_diagram[:, 1, 11] = [6, 1, 7, 0, 5, 2, 4, 3]

        i_bit_mask = np.zeros(3 * bit_length)
        ind        = np.arange(bit_length)
#       order      = dblarr(npoint)
 
        for ip in range(npoint): # check if loop range is valid. 
#            print(ip)
        ## convert to binary
            x_bit_mask = self._btest(x[ip], bit_length - 1, True)
            y_bit_mask = self._btest(y[ip], bit_length - 1, True)
            z_bit_mask = self._btest(z[ip], bit_length - 1, True)

            ## interleave bits
            i_bit_mask[3 * ind + 2] = x_bit_mask
            i_bit_mask[3 * ind + 1] = y_bit_mask
            i_bit_mask[3 * ind + 0] = z_bit_mask

            ## build Hilbert ordering using state diagram
            cstate = 0
            # from bit_length -1 to 0 in descending order.
            for i in range(bit_length -1, -1, -1):
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
                i_bit_mask[3 * i + 2] = self._btest(hdigit, 2, p_all = False)
                i_bit_mask[3 * i + 1] = self._btest(hdigit, 1, p_all = False)
                i_bit_mask[3 * i + 0] = self._btest(hdigit, 0, p_all = False)
                cstate = nstate
                
            ## save Hilbert key as double precision real
            order=[0.0]*npoint
            for i in range(3 * bit_length ):
                b0 = 0
                if (i_bit_mask[i]):
                    b0 = 1
                order[ip] += b0 * (2**i)
#                print(i,ip,i_bit_mask[i], order[ip])
        return order

    def _btest(self, tmp, bit, p_all=False):
        # If p_all = True: it's a vector calculation
        #tmp = i
        nbit = bit
        if (not p_all and tmp != 0): 
            tmp2 = int(np.log2(tmp))+1
            if (tmp2 > nbit): nbit=tmp2
        
        res = [0]*(nbit+1)

        for j in np.arange(nbit +1, 0, -1) - 1:
            res[j] = np.int(tmp / 2**j)
            tmp   -= res[j] * 2**j
        
        if (p_all) :
            return(res)
        else :
            return(res[bit])
#%%

class Info:
    def __init__(self, nout, base):
        self.nout = nout
        self.base = base
        snout = str(self.nout).zfill(5)
        self.fn = base + 'output_' + snout + '/info_' + snout + '.txt'
        self._read_info() # sets ncpu, 

# Weak "internal use" indicator.
# from M import * does not import _blah_blah. 

    def _cal_units(self,arr,rarr):
        # in cgs unit
        kpc     = 3.08e21
        twopi   = 6.2831853e0
        hplanck = 6.6262000e-27
        eV      = 1.6022000e-12
        kB      = 1.38e-16
        clight  = 2.9979250e+10
        Gyr     = 3.1536000e+16
        X       = 0.76
        Y       = 0.24 
        rhoc    = 1.8800000e-29
        mH      = 1.6600000e-24 
        mu_mol  = 1.2195e0
        G       = 6.67e-8
        m_sun   = 1.98892e33
        scale_l = rarr[8]
        scale_d = rarr[9]
        scale_t = rarr[10]
        scale_v = scale_l / scale_t
        scale_T2 = mH/kB * scale_v**2
        scale_nH = X/mH * scale_d
        scale_Z = 1./0.02
        scale_flux = scale_v*scale_d*kpc*kpc*Gyr/m_sun

        self.ncpu = arr[0]
        self.boxtokpc = rarr[0] * scale_l/kpc
        self.tGyr = rarr[1] * scale_t/Gyr
        #self.nout = nout
        self.boxlen = rarr[0]
        self.lmin = arr[2]
        self.lmax = arr[3]
        self.unit_l = scale_l
        self.unit_d = scale_d
        self.unit_t = scale_t
        self.unit_v = scale_v
        self.unit_nH = scale_T2
        self.unit_T2 = scale_Z
        self.kms = scale_v/1e5
        self.unit_flux = scale_d * scale_v * (1e-9*Gyr)*kpc/m_sun
        self.punit_m = scale_d * scale_l**3/2e33
        self.pboxsize = rarr[0] * scale_l/(kpc*1000.)
        self.time = rarr[1]
        self.aexp = rarr[2]
        self.zred = 1/rarr[2]-1
        self.H0 = rarr[3]
        self.om = rarr[4]
        self.ol = rarr[5]
        self.ok = rarr[6]
        self.ob = rarr[7]
        self.msun = scale_d*scale_l**3/2e33

    def keys(self):
        from pprint import pprint
        pprint(vars(self))


    def _read_info(self):
        with open(self.fn) as f:
            arr1 = []
            arr2 = []
#  1 ncpu        =        240
#  2 ndim        =          3
#  3 levelmin    =          7
#  4 levelmax    =         19
#  5 ngridmax    =     800000
#  6 nstep_coarse=      12103
            for i in range(5):
                arr1.append( int(str.split(f.readline(), '=')[1]))
            self.ndim = arr1[1]
            self.nstep_coarse = int(str.split(f.readline(), '=')[1])
            f.readline() # empty line
            for i in range(11):
                arr2.append( float(str.split(f.readline(), '=')[1]) )
            self._cal_units(arr1,arr2)
            self.hilbertkey = np.zeros((2,self.ncpu+1),dtype=np.float32)
            f.readline()
            f.readline()
            f.readline()
            for i in range(self.ncpu):
                keystr = (str.split(f.readline()))[1:]
                self.hilbertkey[:,i]= keystr[0]
                self.hilbertkey[:,i+1]=keystr[1]
#%%      

# AMR header (from IDL)
# 
# ; Initialize header variables
# ncpu_run=0L & ndim=0L & nx=0L & ny=0L & nz=0L
# nlevelmax=0L & ngridmax=0L & nboundary=0L & ngridactual=0L & nstep=0L
# noutput=0L & boxlen=0.0d0 & t=0.0d0
# iout=0L & ifout=0L
# aexp=0.0d0 & hexp=0.0d0 & aexp_old=0.0d0 & epot_tot_int=0.0d0
# epot_tot_old=0.0d0
# const=0.0d0 & mass_tot_0=0.0d0 & rho_tot=0.0d0
# omega_m=0.0d0 & omega_l=0.0d0 & omega_k=0.0d0 & omega_b=0.0d0 & h0=0.0d0
# aexp_ini=0.0d0 & mass_sph=0.0d0
# headf=0L & tailf=0L & numbf=0L & used_mem=0L & used_mem_tot=0L
# 


class Part:
    """ A particle class wih different types of particles, and other meta data.
    DM, star, sink
    """
# There are two types of methods; static method and class method.

# a decorator indicating that the following method does not depend 
# on the actual instance
    @staticmethod 
    def welcome():
        print("It's Part class")
        
    def __init__(self, info, ptypes=None):
# Rule of thumb: initialize everything in __init__ 
# Otherwise, a method may fail to find an object to act on, and so on. 
# This is called 'Object consistency'
        self.ptypes = ptypes
        self.ncpu = 0
        self.nstar = 0
        self.nsink = 0
        self.dm = 0
        self.star = 0
        self.sink = 0
        
        snout = str(info.nout).zfill(5)
        # file name
        self._fnbase = info.base + 'output_' + snout + '/part_' + snout + '.out'
        # header structure
        self._ramses_particle_header = np.dtype([('ncpu', 'i4'), ('ndim', 'i4'), ('npart', 'i4'),
                                   ('randseed', 'i4', (4,)), ('nstar','i4'), 
                                   ('mstar', 'f8'),('mstar_lost', 'f8'), ('nsink', 'i4')])
        self._get_basic_info()
        # Depending on user's choice, generate dm, star, sink classes


#static methods, that don't have access to self

    
    def _get_basic_info(self):
        f = open(self._fnbase + '00001', "rb")
        header = read_header(f, self._ramses_particle_header)
        self.ncpu = header['ncpu']
        self.nstar = header['nstar']
        self.nsink = header['nsink']
        

    def _get_npart_arr(self,cpus):
        print('_get_npart_arr')
        npart_arr = []
        for icpu in cpus:
            with open(self._fnbase + str(icpu).zfill(5), "rb") as f:
                header = read_header(f, self._ramses_particle_header)
                npart_arr.append(header['npart'])

        return npart_arr#, nstar_this, my_mask

    def _get_zoomin(self):
        print(self)
        print(self.has_dm)
        if not (self.has_dm):
            self.add_part("dm")
            self.part.load(self)
            
        # least massive DM particles
        # Or, is there amr level of DM?
        mass = self.part.dm['m']
        print(np.unique(mass))
        min_mass = np.min(mass)

    
    def help():
        print(" Add some helps later on ")

    def load(self,sim,zoom=False,ranges=None):
        # only xyz coordinate is useful with Hilbert space domain decomposition information.
        print("Loading particle...")
        '''
        load(self,zoom=False,xr=[0,1],yr=[0,1],zr=[0,1]):
        self.dm['px'] = np.zeros((self.ndim,self.ndm), dtype='f8')
        self.dm.vel = np.zeros((self.ndim,self.ndm), dtype='f8')
        self.dm.m = np.zeros(self.ndm, dtype='f8')
        # Mass may also be omitted !
        self.dm.id = np.zeros(self.ndm, dtype='f8')
        #self.dm.t = np.zeros(self.ndm, dtype='f8')
        #self.dm.z = np.zeros(self.ndm, dtype='f8') no need to store
        self.star.pos = np.zeros((self.ndim,self.nstar), dtype='f8')
        self.star.vel = np.zeros((self.ndim,self.nstar), dtype='f8')
        self.star.m = np.zeros(self.ndm, dtype='f8')
        self.star.id = np.zeros(self.ndm, dtype='f8')
        self.star.t = np.zeros(self.ndm, dtype='f8')
        self.star.z = np.zeros(self.ndm, dtype='f8')
        '''
        if not(ranges): ranges = sim.ranges
        sim.ranges = ranges        
        
        # Total particle number from selected cpus.
        npart_arr = self._get_npart_arr(sim.cpus)
        print(sim.cpus)
        print('npart_arr:', npart_arr)
        npart_tot = sum(npart_arr)

        # Calculate total number of DM, star, sink from selected cpus.
        # partilce ID, time are needed to distinguish particle type.
        ndm_tot = 0
        nstar_tot = 0
        nsink_tot = 0
        
        for icpu in sim.cpus:
#            if icpu == max(sim.cpus): _end = '\n'
#            else: _end = '\r'
#            print("Loading particles in %d-th cpu output out of %d cpus." % (icpu,len(sim.cpus) ), end=_end)

            f = open(self._fnbase + str(icpu).zfill(5), "rb") #########################  +1 

            header_icpu = read_header(f, self._ramses_particle_header)
            npart_icpu = header_icpu['npart']

            # read position and determine number of particles in the ranges.
            read_fortran(f, np.dtype('f8'), npart_icpu)
            read_fortran(f, np.dtype('f8'), npart_icpu)
            read_fortran(f, np.dtype('f8'), npart_icpu)

            # skip velocity
            read_fortran(f, np.dtype('f8'), npart_icpu)
            read_fortran(f, np.dtype('f8'), npart_icpu)
            read_fortran(f, np.dtype('f8'), npart_icpu)

            # skip mass
            read_fortran(f, np.dtype('f8'), npart_icpu) 

            # read particle id
            id_temp = read_fortran(f, np.dtype('i4'), npart_icpu)             

            # skip refinement
            read_fortran(f, np.dtype('i4'), npart_icpu) 

            # read time
            t_temp = read_fortran(f, np.dtype('f8'), npart_icpu)  
            
            # star particles have positive creation time.
            # mostly positive ID, but young stars 
            # (before SN burst) have negative IDs.
            i_star = (abs(t_temp) > 0.0000001)
            # Dark Matter particles have 0 creatino time, positive ID
            i_dm = np.logical_and(id_temp > 0, t_temp == 0)
            # Sink (BH) particles have negative ID and creation time 0.            
            i_sink = np.logical_and(id_temp < 0, t_temp == 0)
            nstar_icpu = sum(i_star)
            ndm_icpu = sum(i_dm)
            nsink_icpu = sum(i_sink)
            #print(npart_icpu)
            print("icpu:", icpu)            
            print('npart, nStar, nDM, nSink', npart_icpu, nstar_icpu, ndm_icpu, nsink_icpu)
            ndm_tot = ndm_tot + ndm_icpu
            nstar_tot = nstar_tot + nstar_icpu
            nsink_tot = nsink_tot + nsink_icpu
            print("Total")
            print('npart, nStar, nDM, nSink', npart_tot, nstar_tot, ndm_tot, nsink_tot)
            
        # number of darkmatter = npart - nstar - nsink * 2109
        # But!! nstar and nsink is for the whole simulation volume while 
        # npart is for only selected cpus. hmm.        
        
#        self.ndm = sum(npart_arr) - self.nstar - self.nsink * 2109
        
        # Total number of particles stored in the cpus.
        # But particles only within ranges, not in cpus, are eventually returned.
        # So ndm_tot is not going to be the size of DM array.



        #self.ndm = ndm_tot
        #self.nstar = nstar_tot
        #self.nsink = nsink_tot / 2109

        print("Total DM particle %d" % ndm_tot)
        print("Total star particle %d" % nstar_tot)
        print("Total sink particle %d" % nsink_tot)

        # make an array 
        # iterate over files to read in data 
        if (self.has_dm): i_skip_dm = 0
        if (self.has_star): i_skip_star = 0
        if (self.has_sink): i_skip_sink = 0


        for icpu in sim.cpus:
            #if icpu == max(sim.cpus): _end = '\n'
            #else: _end = '\r'
            #print("Loading particles in %d-th cpu output out of %d cpus." % (icpu,len(sim.cpus) ), end=_end)

            f = open(self._fnbase + str(icpu).zfill(5), "rb") #########################  +1 

            header_icpu = read_header(f, self._ramses_particle_header) # skip header
            #header_icpu['nboundary']
            npart_icpu = header_icpu['npart']

            # position
            x_temp = read_fortran(f, np.dtype('f8'), npart_icpu) # row-major 
            y_temp = read_fortran(f, np.dtype('f8'), npart_icpu)
            z_temp = read_fortran(f, np.dtype('f8'), npart_icpu)
            
            if ((ranges[0][0] > 0) or (ranges[0][1] < 1) or
                (ranges[1][0] > 0) or (ranges[1][1] < 1) or
                (ranges[2][0] > 0) or (ranges[2][1] < 1)):

                do_slice=True
                # cut out irrelevant particles
    
                range_ok = np.where(  (ranges[0][0] < x_temp) 
                                    & (ranges[0][1] > x_temp) 
                                    & (ranges[1][0] < y_temp) 
                                    & (ranges[1][1] > y_temp) 
                                    & (ranges[2][0] < z_temp)
                                    & (ranges[2][1] > z_temp) )                
            else: 
                do_slice=False
                range_ok = np.ones_like(x_temp, dtype=np.int)


            # make views to the original array
            px_temp = x_temp[range_ok] 
            py_temp = y_temp[range_ok]
            pz_temp = z_temp[range_ok]
            
            # velocity
            if do_slice:
                #npart_icpu = np.sum(range_ok)
                #npart_icpu = np.product(np.shape(range_ok))
                pass
            
            vx_temp = read_fortran(f, np.dtype('f8'), npart_icpu)[range_ok]
            vy_temp = read_fortran(f, np.dtype('f8'), npart_icpu)[range_ok]
            vz_temp = read_fortran(f, np.dtype('f8'), npart_icpu)[range_ok]

            # mass
            m_temp = read_fortran(f, np.dtype('f8'), npart_icpu)[range_ok]

            # id
            id_temp = read_fortran(f, np.dtype('i4'), npart_icpu)[range_ok]

            # refinement
            ref_temp = read_fortran(f, np.dtype('i4'), npart_icpu)[range_ok] 

            # time
            t_temp = read_fortran(f, np.dtype('f8'), npart_icpu)[range_ok]  

            # metal
            z_temp = read_fortran(f, np.dtype('f8'), npart_icpu)[range_ok] 






        if ( 'dm'  in self.ptypes or  'DM'  in self.ptypes):
            print('')
            print("Loading Dark Matter particles")
            print('')
            self.has_dm = True
            dtype = [('px', '<f8'), ('py', '<f8'), ('pz', '<f8'),
                     ('vx', '<f8'), ('vy', '<f8'), ('vz', '<f8'),
                     ('m', '<f8'), ('id', '<i4'),('ref', '<i4')]
            self.dm = np.recarray(self.ndm, dtype=dtype)#, names = dm_names)
        else: self.has_dm = False

        if ( 'star'  in self.ptypes or  'STAR'  in self.ptypes):
            self.has_star = True
            dtype = [('px', '<f8'), ('py', '<f8'), ('pz', '<f8'),
                     ('vx', '<f8'), ('vy', '<f8'), ('vz', '<f8'),
                     ('m', '<f8'), ('id', '<i4'), ('t', '<f4'), ('z', '<f8')]
            self.star = np.recarray( self.nstar, dtype=dtype)
        else: self.has_star = False

        if ( 'sink'  in self.ptypes or  'SINK'  in self.ptypes):
            self.has_sink = True
            dtype = [('px', '<f8'), ('py', '<f8'), ('pz', '<f8'),
                     ('vx', '<f8'), ('vy', '<f8'), ('vz', '<f8'),
                     ('m', '<f8'), ('id', '<i4')]
            self.sink = np.recarray( self.nsink * 2109, dtype=dtype)
        else: self.has_sink = False




            # distinguish sink / dm / star
            # non star : t == 0
            # sink : t ==0, id < 0

# Copy data to form contiguous arrays of particles.

            if ( 'star'  in self.ptypes or  'STAR'  in self.ptypes):
                i_star = (abs(t_temp) > 0.0000001)
                nstar_icpu = sum(i_star)

                self.star['px'][i_skip_star:i_skip_star + nstar_icpu] = px_temp[i_star]
                self.star['py'][i_skip_star:i_skip_star + nstar_icpu] = py_temp[i_star]
                self.star['pz'][i_skip_star:i_skip_star + nstar_icpu] = pz_temp[i_star]
                self.star['vx'][i_skip_star:i_skip_star + nstar_icpu] = vx_temp[i_star]
                self.star['vy'][i_skip_star:i_skip_star + nstar_icpu] = vy_temp[i_star]
                self.star['vz'][i_skip_star:i_skip_star + nstar_icpu] = vz_temp[i_star]
                self.star['m' ][i_skip_star:i_skip_star + nstar_icpu] = m_temp[i_star]
                self.star['id'][i_skip_star:i_skip_star + nstar_icpu] = id_temp[i_star]
                self.star['t'] [i_skip_star:i_skip_star + nstar_icpu] = t_temp[i_star]
                self.star['z'] [i_skip_star:i_skip_star + nstar_icpu] = z_temp[i_star]
                i_skip_star += nstar_icpu

            #i_dm = id_temp < 0
            i_dm = np.logical_and(id_temp > 0, t_temp == 0)
            i_sink = np.logical_and(id_temp < 0, t_temp == 0)
            ndm_icpu = sum(i_dm)
            if (ndm_icpu != self.ndm):
                print("ndm_icpu:",ndm_icpu)
                print("self.ndm:",self.ndm)
            
            nsink_icpu = sum(i_sink)
            #print('nDM, nSink', ndm_icpu, nsink_icpu)

            #print(ndm_icpu, nsink_icpu)
            # Note that if it's two-division separation,
            # i_dm = t_temp == 0 
            # and then, it's faster to use ~i_dm than to generate another index array.
            if ( 'dm'  in self.ptypes or  'DM'  in self.ptypes):
                self.dm['px'][i_skip_dm:i_skip_dm + ndm_icpu] = px_temp[i_dm]
                self.dm['py'][i_skip_dm:i_skip_dm + ndm_icpu] = py_temp[i_dm]
                self.dm['pz'][i_skip_dm:i_skip_dm + ndm_icpu] = pz_temp[i_dm]
                self.dm['vx'][i_skip_dm:i_skip_dm + ndm_icpu] = vx_temp[i_dm]
                self.dm['vy'][i_skip_dm:i_skip_dm + ndm_icpu] = vy_temp[i_dm]
                self.dm['vz'][i_skip_dm:i_skip_dm + ndm_icpu] = vz_temp[i_dm]
                self.dm['m'] [i_skip_dm:i_skip_dm + ndm_icpu] = m_temp[i_dm]
                self.dm['id'][i_skip_dm:i_skip_dm + ndm_icpu] = id_temp[i_dm]
                self.dm['ref'][i_skip_dm:i_skip_dm + ndm_icpu] = ref_temp[i_dm]
                i_skip_dm += ndm_icpu


            # Which is faster? 
            # i_star[i_dm] as ndm array
            # or i_dm as npart array + i_sink as npart array 
            if ( 'sink'  in self.ptypes or  'SINK'  in self.ptypes):
                self.sink['px'][i_skip_sink:i_skip_sink + nsink_icpu] = px_temp[i_sink]
                self.sink['py'][i_skip_sink:i_skip_sink + nsink_icpu] = py_temp[i_sink]
                self.sink['pz'][i_skip_sink:i_skip_sink + nsink_icpu] = pz_temp[i_sink]
                self.sink['vx'][i_skip_sink:i_skip_sink + nsink_icpu] = vx_temp[i_sink]
                self.sink['vy'][i_skip_sink:i_skip_sink + nsink_icpu] = vy_temp[i_sink]
                self.sink['vz'][i_skip_sink:i_skip_sink + nsink_icpu] = vz_temp[i_sink]
                self.sink['m'] [i_skip_sink:i_skip_sink + nsink_icpu] = m_temp[i_sink]
                self.sink['id'][i_skip_sink:i_skip_sink + nsink_icpu] = id_temp[i_sink]
                i_skip_sink += nsink_icpu


_head_type = np.dtype('i4')

def skip_fortran(f, n=1, verbose=False):

    # A highly efficient way of reading binary data with a known data-type
    alen = np.fromfile(f, _head_type, 1) # == skip
    if verbose:
        print("FORTRAN block length %d!=%d" % (alen))
        
    mod_check = alen % 4
    if mod_check != 0:
        print("Array size is not a multiple of 4")
    
    n = int(alen/4)

    data = np.fromfile(f, _head_type, n) 
    #print('check',data)
    alen = np.fromfile(f, _head_type, 1)


def read_fortran(f, dtype, n=1):
    if not isinstance(dtype, np.dtype):
        dtype = np.dtype(dtype)

    length = n * dtype.itemsize
    
    alen = np.fromfile(f, _head_type, 1) # == skip
    
    if alen != length:
        raise IOError("Unexpected FORTRAN block length %d!=%d" % (alen, length))

    
    data = np.fromfile(f, dtype, n) # Actual data
    
    alen = np.fromfile(f, _head_type, 1)
    if alen != length:
        raise IOError("Unexpected FORTRAN block length (tail) %d!=%d" % (alen, length))

    return data

def read_header(f, dtype):
    q = np.empty(1, dtype=dtype)
    for i in range(len(dtype.fields)):
        data = read_fortran(f, dtype[i] )      
    
        # I really don't understand why the following acrobatic should
        # be necessary, but q[0][i] = data[0] doesn't copy arrays properly
        if np.issubdtype(dtype[i], np.string_):
            q[0][i] = data
        elif hasattr(data[0], "__len__"):
            q[0][i][:] = data[0]
        else:
            q[0][i] = data[0]
        # if string, return the whole array

    return q[0]


def read_header_string(f, dtype):
    alen = np.fromfile(f, _head_type, 1) # == skip
    length = dtype.itemsize

    if alen != length:
        raise IOError("Unexpected FORTRAN block length %d!=%d" % (alen, length))

    data = np.fromfile(f, dtype, 1) # Actual data
    #print('%s',self.ordering)
    alen = np.fromfile(f, _head_type, 1)
    if alen != length:
        raise IOError("Unexpected FORTRAN block length (tail) %d!=%d" % (alen, length))
    return data

#%%

# How about reading stars and dms per cpu output? 
# No way to know the number of stars in individual output beforehand. 
# Read all files once more only to get the size of stars? Maybe no!
# Or, append numpy array hundreds or thousands of times to construct one big array? Well...
# 
# For the moment, I'll get the total number of stars and non-stars first. 
# Then, allocate momeries of two populations.
# 
# Question.
# Which of the three is the best for performance's sake? mask array, view, or index array?
# 
# 
# + 
# Later, add Hilber-domain decomposition consideration. - done


class AmrHeader():
    def __init__(self):
        pass
    
    def _read_amr_header(self,f):
    # Make this visible from outside, and more general
    # this can be used everytime you need to skip header
    # where value assigning is unnecessary
    # recieve file object f rather than opening one internally.
    #
    # Are all header entries global quantaties? 
    # or varies with cpus?
    # Local values must be separated from global values.
                

        # Global
        h1 = read_header(f, np.dtype([('ncpu', 'i4'), ('ndim', 'i4'), ('ng', 'i4', (3,)),
                                       ('nlevelmax', 'i4'), ('ngridmax','i4'), ('nboundary', 'i4'),
                                       ('ngrid', 'i4'), ('boxlen', 'f8'), ('outputs', 'i4', (3,))]))
        
        # Basic information that are required to read the header further.
        ncpu = h1['ncpu']
        nnouts = h1['outputs'][0]
        nlevelmax = h1['nlevelmax']
        ncoarse = np.product(h1['ng'][:]) # multiply all elements in an array        
        
        # Global
        h2 = read_header(f, np.dtype([('tout', 'f8', (nnouts,)), ('aout', 'f8', (nnouts,)), ('t', 'f8'),
                                      ('dtold', 'f8', (nlevelmax,)), ('dtnew', 'f8', (nlevelmax,)),
                                      ('nsteps', 'i4', (2,)), 
                                      ('cmr','f8', (3,)),
                                      ('omlkbhal','f8', (7,)),
                                      ('expepot','f8', (5,)),
                                      ('mass_sph', 'f8'),
                                      ('headl', 'i4', (nlevelmax,ncpu,)),
                                      ('taill', 'i4',(nlevelmax,ncpu,)),
                                      ('numbl', 'i4',(nlevelmax,ncpu,))]))
        skip_fortran(f)
    # headl
    #
    
    # When reading 2D array, beware that fortran file is written in 
    # ???-major order but python will save it in ???-major order
    
    # For example, 
    # numbl is (ncpu x nlevelmax) array in fortran
    # and is accessed by numbl[icpu,ilevel]
    
    # But in Python, this should be
    # 
    
    # where is the 170 coming from?                             
        if (h1['nboundary'] > 0):
            h3 = read_header(f, np.dtype([('headb', 'i4',(nlevelmax,ncpu,)),
                                          ('tailb', 'i4',(nlevelmax,ncpu,)),
                                          ('numbb', 'i4',(nlevelmax,ncpu,))])) 
            self.headb = h3['headb']
            self.tailb = h3['tailb']
            self.numbb = h3['numbb']

        h4 = read_header(f, np.dtype([('htnm1m2', 'i4', (5,)),
                                      ('ordering', 'a128'),
                                      ('bound_key', 'i8', (ncpu+1,)),
                                      ('son', 'i4', (ncoarse,)),
                                      ('flag1', 'i4', (ncoarse,)),
                                      ('cpu_map', 'i4', (ncoarse,))]))

        # if assign 
        self.ncpu = h1['ncpu']
        self.ndim = h1['ndim']
        self.ng = h1['ng']
        self.nlevelmax = h1['nlevelmax']
        self.ngridmax = h1['ngridmax']
        self.nboundary = h1['nboundary']
        self.ngrid = h1['ngrid']
        self.boxlen = h1['boxlen']
        self.nnouts = h1['outputs'][0]
        self.iout = h1['outputs'][1]
        self.ifout = h1['outputs'][2]
        
        self.tout = h2['tout']
        self.aout = h2['aout']
        self.t = h2['t']
        self.dtold = h2['dtold']
        self.dtnew = h2['dtnew']       
        self.const = h2['cmr'][0]
        self.mass_tot0 = h2['cmr'][1]
        self.rho_tot = h2['cmr'][2]
        self.nstep = h2['nsteps'][0] 
        self.nstep_coarse = h2['nsteps'][1]
        
        self.Om       = h2['omlkbhal'][0]
        self.Ol       = h2['omlkbhal'][1]
        self.Ok       = h2['omlkbhal'][2]
        self.Ob       = h2['omlkbhal'][3]
        self.h0       = h2['omlkbhal'][4]
        self.aexp_ini = h2['omlkbhal'][5]
        self.boxlen   = h2['omlkbhal'][6]
        
        self.aexp         = h2['expepot'][0]
        self.hexp         = h2['expepot'][1]
        self.aexp_old     = h2['expepot'][2]
        self.epot_tot_ini = h2['expepot'][3]
        self.epot_tot_old = h2['expepot'][4]
        
        self.mass_sph = h2['mass_sph']
        
        self.headl  = h2['headl']
        self.taill  = h2['taill']
        self.numbl  = h2['numbl']
#        self.numbot = h2['numbot'] # This value has been skipped
        
        self.headf        = h4['htnm1m2'][0]
        self.tailf        = h4['htnm1m2'][1]
        self.numbf        = h4['htnm1m2'][2]
        self.used_mem     = h4['htnm1m2'][3]
        self.used_mem_tot = h4['htnm1m2'][4]
        
        self.ordering  = h4['ordering']
        self.bound_key = h4['bound_key']
        self.son       = h4['son']
        self.flag1     = h4['flag1']
        self.cpu_map   = h4['cpu_map']


class Amr():

    def __init__(self, info):
        snout = str(info.nout).zfill(5)
        # file name
        self._fnbase = info.base + 'output_' + snout + '/amr_' + snout + '.out'
        f = open(self._fnbase + '00001', "rb") # 00001 = first of the cpu list
        
        # header structure
        self.header = AmrHeader()
        self.header._read_amr_header(f)
        f.close()

    def keys(self):
        from pprint import pprint
        pprint(vars(self))
       
   
    def _load_mesh(f, ndim = 3):
        for i in np.arange(ndim):
            read_fortran(f, np.dtype('f8'))
    
    def load(self,sim,verbose=False):
        
# need some information from info file.

# The building block of FTT AMR is an oct, which is a group of 8 cells and data
# either associate with cells or the oct(called 'mesh'8 in IDL analysis routine).
# An Oct consists of (level, xyz coordinates, poiner to the parent,
# pointer to 6 neighbouring parents, pointer to child oct)

# In RAMSES, cpu map and refinement map is additionaly needed to restart a simulation.
# In RAMSES, octs of the same level are written at once. (So you need to loop over ilevel)
# 

#####
##### Think about when you will need to load the amr file. 
##### It's rather unclear! 

        # global header variables are available.
        icpu = 0 ## icpu 
        cpus = sim.cpus
        ndim = self.header.ndim
        ngridtot = 0
        
        ncpu = self.header.ncpu
        nboundary = self.header.nboundary
        nlevelmax = self.header.nlevelmax
        
        # Size of arrays
        
        if icpu == 0:
            listmax = ncpu 
        elif icpu > 0:
            listmax = ncpu + nboundary
                
        ngrid = np.zeros((nlevelmax,listmax), dtype=np.int32)
        
        llist = 0
        for jcpu in cpus:
#            if jcpu == max(cpus): _end = '\n'
#            else: _end = '\r'
#            print("Loading AMR grids in %d-th cpu output out of %d cpus." % (jcpu,max(cpus)), end=_end)
            print(self._fnbase + str(jcpu).zfill(5))
            f = open(self._fnbase + str(jcpu).zfill(5), "rb") #########################  +1 

            # read header 
            header_icpu = AmrHeader()
            header_icpu._read_amr_header(f)
            
            numbl = header_icpu.numbl
            if nboundary > 0 : 
                numbb = header_icpu.numbb # None if nboundary = 0
            else:
                numbb = np.zeros(np.shape(numbl))
            # ngrid != header_icpu.ngrid 
            #            
            
            ngridtot = 0
            
            kcpumin = 1
            kcpumax = nboundary + ncpu
            nlevel = 0        
            
            for ilevel in np.arange(nlevelmax):
                for kcpu in np.arange(kcpumin, kcpumax + 1):

                    if (kcpu <= ncpu): 
                        ng = numbl[ilevel][kcpu-1]
                    else: 
                        ng = numbb[ilevel][kcpu - ncpu -1]
                    if icpu == 0:
                        if kcpu == jcpu: 
                            ngrid[ilevel][kcpu-1] = ng
                    else: 
                        ngrid[ilevel][kcpu-1] = ng 
                    
                    if (ng > 0):
                        ngridtot = ngridtot + ng
                        if (verbose): 
                            print("Level %2d has %6d grids in proc %4d" % (ilevel+1,ng,kcpu))
                        print(f.tell())
                        nlevel = nlevel + 1 # number of valid (ng > 0) levels
                        

                        # Read actual data
                        i_skip = 0
                        ind_current = read_fortran(f, np.dtype('i4'), ng) # row-major 
                        ind_next = read_fortran(f, np.dtype('i4'), ng)
                        ind_prev = read_fortran(f, np.dtype('i4'), ng)
                        for idim in range(ndim): # gird center
                            xx = read_fortran(f, np.dtype('f8'), ng)
                        # father index
                        read_fortran(f, np.dtype('i4'), ng)
                        
                        for idim in range(2*ndim): # neighbour index
                            read_fortran(f, np.dtype('i4'), ng)
                            
                        for idim in range(2**ndim): # son index
                            read_fortran(f, np.dtype('i4'), ng)
                            
                        for idim in range(2**ndim): # cpu map
                            read_fortran(f, np.dtype('i4'), ng)
                            
                        for idim in range(2**ndim): # ref map
                            read_fortran(f, np.dtype('i4'), ng)

                        if icpu == 0 :
                            if (kcpu == jcpu):
                                pc = 1
                            
                        llist += 1
                        
            #self.amr['xx'][i_skip:i_skip + ndm_icpu] = pos_temp[0,i_dm]
            f.close
            
## So, What will you return? 
#%%

class Hydro:
    class HydroHeader():
        def __init__(self):
            pass
    
    def __init__(self, info):
        # self = Hydro
        snout = str(info.nout).zfill(5)
        # file name
        self._fnbase = info.base + 'output_' + snout + '/hydro_' + snout + '.out'       
        self._get_basic_info()
    
    def _get_basic_info(self):
        f = open(self._fnbase + '00001', "rb")
        self.header = self.HydroHeader()
        self._read_hydro_header(f)
    
    def _read_hydro_header(self,f):
        # Global
        h1 = read_header(f, np.dtype([('ncpu', 'i4'), ('nvarh', 'i4'), ('ndim', 'i4'),
                                       ('nlevelmax', 'i4'), ('nboundary', 'i4'),
                                       ('gamma','f8')]))
        print('ncpu',h1['ncpu'])
        print('nvarh',h1['nvarh'])
        print('gamma',h1['gamma'])
        if h1['nboundary'] == 0: print(' Periodic boundary condition')

        # if assign 
        self.header.ncpu = h1['ncpu']
        self.header.ndim = h1['ndim']
        self.header.nlevelmax = h1['nlevelmax']
        self.header.nboundary = h1['nboundary']
        self.header.gamma = h1['gamma']
        self.header.nvarh = h1['nvarh']

            
    def amr2cell(self, sim, lmax = None):
        # if both amr and hydro does not exist, create. 
        
       
        #################
        verbose = True
        icpu = 0 ## icpu 
        cpus = sim.cpus
        ndim = sim.info.ndim
        print(ndim)
        twotondim = 2**ndim
        ngridtot = 0
        
        nvarh = self.header.nvarh
        
        ncpu = sim.amr.header.ncpu
        nboundary = sim.amr.header.nboundary
        nlevelmax = sim.amr.header.nlevelmax
        if lmax == None: lmax = nlevelmax
            
        xmi = sim.ranges[0][0]
        xma = sim.ranges[0][1]
        ymi = sim.ranges[1][0]
        yma = sim.ranges[1][1]
        zmi = sim.ranges[2][0]
        zma = sim.ranges[2][1]

        print( ' >>> working resolution (lmax) =', lmax)
        
        # Size of arrays
        
        if icpu == 0:
            listmax = ncpu 
        elif icpu > 0:
            listmax = ncpu + nboundary
                
        ngrid = np.zeros((nlevelmax,listmax), dtype=np.int32)
        nx = sim.amr.header.ng[0]
        ny = sim.amr.header.ng[1]
        nz = sim.amr.header.ng[2]
        xbound=[nx/2,ny/2,nz/2]
#        llist = 0
        nlinemax = 20000000
        for icpu in cpus:
            famr = open(sim.amr._fnbase + str(icpu).zfill(5), "rb") 

            # read header 
            a_header_icpu = AmrHeader()
            a_header_icpu._read_amr_header(famr)

            
            ngrid[0:nlevelmax][0:ncpu]=a_header_icpu.numbl 
            # because hydro file is written in fortran order 
            # looping over 
            #
            # for ilevel
            #   for icpu
            
            
            
    

            fhydro = open(sim.hydro._fnbase + str(icpu).zfill(5), "rb") 
            self._read_hydro_header(fhydro) 

            # Initialize pointer array
            
            xt = np.zeros(nlinemax, dtype=np.float64)
            yt = np.zeros(nlinemax, dtype=np.float64)
            zt = np.zeros(nlinemax, dtype=np.float64)
            dxt = np.zeros(nlinemax, dtype=np.float64)
            vart = np.zeros((nvarh,nlinemax), dtype=np.float64)
            
            xc = np.zeros((8,3), dtype=np.float64)
            
            for ilevel in np.arange(lmax):
                # geometry
                dx  = 0.5**ilevel # grid size of the current level.
                dx2 = dx/2 # half of the grid size of the current level.
                nx_full = 2**ilevel # maximum possible number of grid cells at the given level.
                ny_full = 2**ilevel
                nz_full = 2**ilevel
    # division results in float!
                for ind in np.arange(1,twotondim): #  1 to 8
                    iz = int((ind-1)/4 )
                    iy = int((ind-1-4*iz)/2)
                    ix = int((ind-1-2*iy-4*iz))
                    xc[ind-1][0]=(ix-0.5)*dx # grid center
                    xc[ind-1][1]=(iy-0.5)*dx
                    xc[ind-1][2]=(iz-0.5)*dx

                # Allocate work arrays
                ngrida = ngrid[ilevel-1][icpu-1]   #   
                    
    #           grid[ilevel].ngrid=ngrida
                if ngrida > 0 :
                    # note that index stars from 0 to reduce memory use
                    xg  = np.zeros((ndim, ngrida), dtype=np.float64)
                    son = np.zeros((twotondim, ngrida), dtype=np.int32)
                    var = np.zeros((nvarh, twotondim, ngrida,), dtype=np.float64)  
                    x   = np.zeros((ndim, ngrida), dtype=np.float64) 
                    rho = np.zeros(ngrida, dtype=np.float64)
                     
                for jcpu in range(1, nboundary+ncpu):
#                    print(icpu,ilevel,kcpu,ngrid[ilevel-1][icpu-1],ngrid[ilevel-1][kcpu-1])
                    if ngrid[ilevel-1,jcpu-1] == 0 :
                        tmp1 = read_fortran(fhydro, np.dtype('i4') )
                        tmp2 = read_fortran(fhydro, np.dtype('i4') )
                        if (verbose): 
                            print(ngrid[ilevel-1,jcpu-1],"Level %2d has %6d grids in proc %4d" % (tmp1,tmp2,jcpu))
#                        pass
                    else :
                        ng = ngrid[ilevel-1][jcpu-1]
                        if (verbose): 
                            print("Level %2d has %6d grids in proc %4d" % (ilevel+1,ng,jcpu))
                        print(famr.tell())
                        read_fortran(famr, np.dtype('i4'), ng) 
                        print(famr.tell())
                        read_fortran(famr, np.dtype('i4'), ng)
                        print(famr.tell())
                        read_fortran(famr, np.dtype('i4'), ng)
                        print(famr.tell())

                        for idim in range(ndim): # gird center
                            temp = read_fortran(famr, np.dtype('f8'), ng)
                            if jcpu == icpu:
                                xg[:][idim] = temp
                            print(famr.tell())


                        #skip
                        read_fortran(famr, np.dtype('i4'), ng) # father index
                        print(famr.tell())
                        #skip
                        for idim in range(2*ndim): # neighbour index
                            read_fortran(famr, np.dtype('i4'), ng)

                        for idim in range(2**ndim): # son index
                            temp = read_fortran(famr, np.dtype('i4'), ng)
                            if jcpu == icpu:
                                son[:][idim] = temp

                        # skip
                        for idim in range(2**ndim):
                            read_fortran(famr, np.dtype('i4'), ng) # cpu map
                        # skip
                        for idim in range(2**ndim):
                            read_fortran(famr, np.dtype('i4'), ng)  # ref map

                        # Hydro


                        for ind in range(twotondim) :
                            for ivar in range(nvarh) :
                                print(ivar)
                                temp = read_fortran(fhydro, np.dtype('f8'), ng) # aaa
                                print(fhydro.tell())
                                if jcpu == icpu:
                                    var[:][ind][ivar] = temp
                                     ##### fix the order. this won't be efficient.
                    
                    if ngrida == 0 :
                        pass
                    else :
                        print(np.shape(x),np.shape(xg),np.shape(xc),np.shape(xbound))
                        
                        
                        for ind in range(1,twotondim):
                            x[0][:] = xg[0][:] + xc[ind-1][0]-xbound[0]
                            x[1][:] = xg[0][:] + xc[ind-1][1]-xbound[1]
                            x[2][:] = xg[0][:] + xc[ind-1][2]-xbound[2]

                            ref = np.zeros(ngrida, dtype=np.int16)
                            if ilevel < lmax:
                                val = np.where(son[ind -1][:] > 0)
                                if sum(val) > 0 : ref[val] = 1

                            ok_cell = np.where( ref == 0  and 
                                    x[0,:] + dx2 >= xmi and
                                    x[1,:] + dx2 >= ymi and
                                    x[2,:] + dx2 >= zmi and
                                    x[0,:] - dx2 <= xma and
                                    x[1,:] - dx2 <= yma and
                                    x[2,:] - dx2 <= zma )
                            Nok_cell = sum(ok_cell)

                        if Nok_cell > 0 :
                            nbegin = nline
                            nend = nline + Nok_cell - 1

                            # If existing array size if not enough. 
                            # Any workaround?
                            if nend > nlinemax :
                                xt   = np.append(xt, np.zeros(nlinemax, dtype=np.float64))                                             
                                yt   = np.append(yt, np.zeros(nlinemax, dtype=np.float64))
                                zt   = np.append(zt, np.zeros(nlinemax, dtype=np.float64))
                                dxt   = np.append(dxt, np.zeros(nlinemax, dtype=np.float64))
                                vart = np.append(vart, np.zeros((nvarh,nlinemax), dtype=np.float64))
    #                            vart = np.zeros((nlinemax*2,nvarh), dtype=np.float64)
    #                            
    #                            for jj in range(nvarh-1) : 
    #                                vart[0:nlinemax-1,jj] = vart2[:,jj]
    #                            nlinemax*=2
    #                            vart2=0
                            xt[nbegin:nend]    = x[ok_cell,0]
                            yt[nbegin:nend]    = x[ok_cell,1]
                            zt[nbegin:nend]    = x[ok_cell,2]
                            dxt[nbegin:nend]   = dx[ok_cell]
                            vart[nbegin:nend,:]= reform[var[ok_cell,ind-1,:]]

                            nline = nend+1


            # loop over ilevel
            famr.close()
            fhydro.close()
                                
