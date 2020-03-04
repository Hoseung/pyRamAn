import numpy as np
import struct 

def get_binary_data(fmt="",
                    ninteg=0,
                    nlines=0,
                    nfloat=0,
                    nstrin=0,
                    nquadr=0,
                    nlongi=0,
                    content=None,
                    correction=0):

    offset = 4*ninteg + 8*(nlines+nfloat+nlongi) + nstrin + nquadr*16 + 4 + correction
    byte_size = {"i":4,"d":8,"q":8}
    if len(fmt) == 1:
        mult = 1
    else:
        mult = eval(fmt[0:len(fmt)-1])
    pack_size = mult*byte_size[fmt[-1]]

    return struct.unpack(fmt, content[offset:offset+pack_size])

def amr2cell_py(self, lmax=None,
                  icpu=0, cpu=False,
                  read_level = False,
                  read_ref = False,
                  verbose=False, debug=False, nvar_read=6):

      read_cpu = cpu
      # if both amr and hydro does not exist, create.
      #################
      cpus = self.cpus
      ndim = self.info.ndim
      if not verbose: print("CPU list:", cpus)
      twotondim = 2**ndim
      nvarh = self.header.nvarh
      ncpu = self.header.ncpu
      nboundary = self.header.nboundary
      nlevelmax = self.header.nlevelmax
      if lmax is None:
          lmax = nlevelmax

      print(' >>> working resolution (lmax) =', lmax)

      # Set ranges
      xmi = self.ranges[0][0]
      xma = self.ranges[0][1]
      ymi = self.ranges[1][0]
      yma = self.ranges[1][1]
      zmi = self.ranges[2][0]
      zma = self.ranges[2][1]

      # Size of arrays
      if icpu == 0:
          listmax = ncpu
      elif icpu > 0:
          listmax = ncpu + nboundary
      ngridarr = np.zeros((nlevelmax, listmax), dtype=np.int32)

      nx = self.amr.header.ng[0]
      ny = self.amr.header.ng[1]
      nz = self.amr.header.ng[2]
      xbound = [int(nx/2), int(ny/2), int(nz/2)]

      # Allocate work arrays
      twotondim = 2**self.info.ndim
      xcent = np.zeros([8,3],dtype=np.float64)
      xg    = np.zeros([self.info.ngridmax,3],dtype=np.float64)
      son   = np.zeros([self.info.ngridmax,twotondim],dtype=np.int32)
      var   = np.zeros([self.info.ngridmax,twotondim,nvar_read+6],dtype=np.float64)
      xyz   = np.zeros([self.info.ngridmax,twotondim,self.info.ndim],dtype=np.float64)
      ref   = np.zeros([self.info.ngridmax,twotondim],dtype=np.bool)

      iprog = 1
      istep = 10
      ncells_tot = 0

      # Control which variable to read here
      var_read = [True]*nvar_read

      work_dir = self.info.base + '/' + self.out_dir + 'output_{:05d}/'.format(self.info.nout)

      data_pieces = dict()
      npieces = 0

      read_header = True

      ncoarse = np.product(self.amr.header.ng)
      nboundary = self.amr.header.nboundary
      key_size = self.amr.header.key_size
      noutput = self.amr.header.nnouts

      ngridlevel = np.zeros([self.info.ncpu_tot+nboundary,self.info.lmax],dtype=np.int32)

      for kcpu in range(1, self.info.ncpu_tot+1):
          # Read binary AMR file
          fn_amr = work_dir +"amr_{:05d}.out{:05d}".format(self.info.nout, kcpu)
          if kcpu in cpus or kcpu == 1:
              with open(fn_amr, mode='rb') as famr:
                  amrContent = famr.read()

                  # Read binary HYDRO file
                  fn_hydro = fn_amr.replace("amr", "hydro")
                  with open(fn_hydro, mode='rb') as fhydro: # b is important -> binary
                      hydroContent = fhydro.read()
                      # All hydro files have the same header, right?
                      if read_header:
                          fhydro.seek(0)
                          self._read_hydro_header(fhydro)
                          read_header = False

          #print("Reading....... ", kcpu)

          nstrin = 0
          #ngridlevel = header_icpu.numbl

          # Read the number of grids
          ninteg = 14+(2*self.info.ncpu_tot*self.info.lmax)
          nfloat = 18+(2*noutput)+(2*self.info.lmax)
          nlines = 21
          ngridlevel[:self.info.ncpu_tot,:] = np.asarray(get_binary_data(fmt="%ii"%(self.info.ncpu_tot*self.info.lmax),\
           content=amrContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat)).reshape(self.info.lmax,self.info.ncpu_tot).T


          # Offset for AMR
          ninteg1 = 14+(3*self.info.ncpu_tot*self.info.lmax)+(10*self.info.lmax)+(3*nboundary*self.info.lmax)+5+3*ncoarse
          nfloat1 = 18+(2*noutput)+(2*self.info.lmax)
          nlines1 = 21+2+3*min(1,nboundary)+1+1+1+3
          nstrin1 = 128 + key_size

          # Offset for HYDRO
          ninteg2 = 5
          nfloat2 = 1
          nlines2 = 6
          nstrin2 = 0

          # Loop over levels
          for ilevel in range(lmax):

              # Geometry
              dxcell=0.5**(ilevel+1)
              dx2=0.5*dxcell
              for ind in range(twotondim):
                  iz=int((ind)/4)
                  iy=int((ind-4*iz)/2)
                  ix=int((ind-2*iy-4*iz))
                  xcent[ind,0]=(float(ix)-0.5)*dxcell
                  xcent[ind,1]=(float(iy)-0.5)*dxcell
                  xcent[ind,2]=(float(iz)-0.5)*dxcell

              # Cumulative offsets in AMR file
              ninteg_amr = ninteg1
              nfloat_amr = nfloat1
              nlines_amr = nlines1
              nstrin_amr = nstrin1

              # Cumulative offsets in HYDRO file
              ninteg_hydro = ninteg2
              nfloat_hydro = nfloat2
              nlines_hydro = nlines2
              nstrin_hydro = nstrin2

              # Loop over domains
              for jcpu in range(nboundary+self.info.ncpu_tot):

                  ncache = ngridlevel[jcpu, ilevel]#header_icpu.numbl[ilevel,jcpu]

                  # Skip two lines of integers
                  nlines_hydro += 2
                  ninteg_hydro += 2

                  if ncache > 0:

                      if jcpu == kcpu-1:
                          # xg: grid coordinates
                          ninteg = ninteg_amr + ncache*3
                          nfloat = nfloat_amr
                          nlines = nlines_amr + 3
                          nstrin = nstrin_amr
                          for n in range(self.info.ndim):
                              offset = 4*ninteg + 8*(nlines+nfloat+n*(ncache+1)) + nstrin + 4
                              xg[:ncache,n] = struct.unpack("%id"%(ncache),           amrContent[offset:offset+8*ncache])

                          # son indices
                          ninteg = ninteg_amr + ncache*(4 + 2 * self.info.ndim)
                          nfloat = nfloat_amr + ncache*self.info.ndim
                          nlines = nlines_amr + 4 + 3*self.info.ndim
                          nstrin = nstrin_amr
                          for ind in range(twotondim):
                              offset = 4*(ninteg+ind*ncache) + 8*(nlines+nfloat+ind) + nstrin + 4
                              son[:ncache,ind] = struct.unpack("%ii"%(ncache), amrContent[offset:offset+4*ncache])
                              # var: hydro variables
                              jvar = 0
                              for ivar in range(nvarh):
                                  if var_read[ivar]:
                                      offset = 4*ninteg_hydro + 8*(nlines_hydro+nfloat_hydro+(ind*nvarh+ivar)*(ncache+1)) + nstrin_hydro + 4
                                      var[:ncache,ind,jvar] = struct.unpack("%id"%(ncache), hydroContent[offset:offset+8*ncache])
                                      jvar += 1

                              # var: coordinates and cell sizes
                              var[:ncache,ind,-6] = float(ilevel+1) # level
                              for n in range(self.info.ndim):
                                  xyz[:ncache,ind,n] = xg[:ncache,n] + xcent[ind,n]-xbound[n]
                                  var[:ncache,ind,-5+n] = xyz[:ncache,ind,n]*self.info.boxlen
                              var[:ncache,ind,-2] = dxcell*self.info.boxlen
                              var[:ncache,ind,-1] = kcpu # CPU number
                              # ref: True if the cell is unrefined
                              ref[:ncache,ind] = np.logical_not(np.logical_and(son[:ncache,ind] > 0,ilevel < lmax-1))

                          # Select only the unrefined cells that are in the region of interest
                          if self.info.ndim == 1:
                              cube = np.where(np.logical_and(ref[:ncache,:], \
                                              np.logical_and((xyz[:ncache,:,0]+dx2)>=xmi, \
                                                             (xyz[:ncache,:,0]-dx2)<=xma)))
                          elif self.info.ndim == 2:
                              cube = np.where(np.logical_and(ref[:ncache,:], \
                                              np.logical_and((xyz[:ncache,:,0]+dx2)>=xmi, \
                                              np.logical_and((xyz[:ncache,:,1]+dx2)>=ymi, \
                                              np.logical_and((xyz[:ncache,:,0]-dx2)<=xma, \
                                                             (xyz[:ncache,:,1]-dx2)<=yma)))))
                          elif self.info.ndim == 3:
                              cube = np.where(np.logical_and(ref[:ncache,:], \
                                              np.logical_and((xyz[:ncache,:,0]+dx2)>=xmi, \
                                              np.logical_and((xyz[:ncache,:,1]+dx2)>=ymi, \
                                              np.logical_and((xyz[:ncache,:,2]+dx2)>=zmi, \
                                              np.logical_and((xyz[:ncache,:,0]-dx2)<=xma, \
                                              np.logical_and((xyz[:ncache,:,1]-dx2)<=yma, \
                                                             (xyz[:ncache,:,2]-dx2)<=zma)))))))
                          else:
                              print("Bad number of dimensions")
                              return 0

                          cells = var[cube]
                          ncells = np.shape(cells)[0]
                          if ncells > 0:
                              ncells_tot += ncells
                              npieces += 1
                              # Add the cells in the master dictionary
                              data_pieces["piece"+str(npieces)] = cells

                      # Now increment the offsets while looping through the domains
                      ninteg_amr += ncache*(4+3*twotondim+2*self.info.ndim)
                      nfloat_amr += ncache*self.info.ndim
                      nlines_amr += 4 + 3*twotondim + 3*self.info.ndim


                      nfloat_hydro += ncache*twotondim*nvarh
                      nlines_hydro += twotondim*nvarh

              # Now increment the offsets while looping through the levels
              ninteg1 = ninteg_amr
              nfloat1 = nfloat_amr
              nlines1 = nlines_amr
              nstrin1 = nstrin_amr

              ninteg2 = ninteg_hydro
              nfloat2 = nfloat_hydro
              nlines2 = nlines_hydro
              nstrin2 = nstrin_hydro

      all_array = np.concatenate(list(data_pieces.values()), axis=0)

      dtype_cell = {'pos': (('<f8', (3,)), 0),
                      'x': (('<f8', 1), 0),
                      'y': (('<f8', 1), 8),
                      'z': (('<f8', 1), 16),
                     'dx': (('<f8', 1), 24),
                   'var0': (('<f8', 1), 32),
                    'rho': (('<f8', 1), 32),
                    'vel': (('<f8', (3,)), 40),
                     'vx': (('<f8', 1), 40),
                     'vy': (('<f8', 1), 48),
                     'vz': (('<f8', 1), 56),
                   'var1': (('<f8', 1), 40),
                   'var2': (('<f8', 1), 48),
                   'var3': (('<f8', 1), 56),
                   'var4': (('<f8', 1), 64),
                   'temp': (('<f8', 1), 64),
                   'var5': (('<f8', 1), 72),
                  'metal': (('<f8', 1), 72)}
      dt_off = 72
      if read_cpu:
          dtype_cell.update({'cpu': (('<f8',1),dt_off+8)})
          dt_off += 8
      if read_level:
          dtype_cell.update({'level': (('<f8',1),dt_off+8)})
          dt_off +=8
      if read_ref:
          dtype_cell.update({'ref': (('bool',1),dt_off+8)})


      # This part
      self.cell = np.zeros(len(all_array), dtype=dtype_cell)

      for i in range(nvarh):
          self.cell['var' + str(i)] = all_array[:,i]

      if read_level:
          self.cell['level'] = all_array[:,nvar_read+0]

      if read_cpu:
          self.cell['cpu'] = all_array[:,nvar_read+5]

      self.cell['x'] = all_array[:,nvar_read+1]
      self.cell['y'] = all_array[:,nvar_read+2]
      self.cell['z'] = all_array[:,nvar_read+3]
      self.cell['dx'] = all_array[:,nvar_read+4]

      if read_ref:
          self.cell["ref"] = all_array[:,7]

      #return master_data_array
