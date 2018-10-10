# Copyright (c) 2018--, Aura Obreja and the GalacticStructureFinder (gsf) contributors.
# GalacticStructureFinder decomposes simulated galaxies based on their stellar kinematics.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# This program uses third-party libraries:
# - numpy (http://www.numpy.org/)
# - scipy (https://www.scipy.org/)
# - matplotlib (https://matplotlib.org/)
# - pynbody (https://github.com/pynbody/pynbody)
# - scikit-learn (http://scikit-learn.org/stable/)
# Please refer to the original licensing conditions for each of these third party software libraries.
#
# If you use this program or parts of it, please cite the original article presenting gsf:
# Obreja, Maccio, Moster et al MNRAS 2018, "Introducing galactic structure finder: the
# multiple stellar kinematic structures of a Milky Way mass galaxy"
#

import pickle, os, glob, time
from sys import stdin
import random
import numpy as np
import scipy
import scipy.interpolate
import pynbody
import pynbody.filt as filt
import pynbody.units as units
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib.colors import LogNorm
from sklearn.preprocessing import scale
from sklearn.mixture import GMM
import twobody

grav_const=4.302e-6
plt.switch_backend('agg')


def secondsToStr(t):
    return "%d:%02d:%02d.%03d" % \
        reduce(lambda ll,b : divmod(ll[0],b) + ll[1:],
            [(t*1000,),1000,60,60])

def wgsf(snaphot_file, out_dir='./', number_of_clusters=2, covariance_type='full', whiten_data=None,
                halo_id=1, radius_align=0.1, align_with='baryon', n_init=1, plot=True, verbose=True):

    arp = snaphot_file.split('/')
    if len(arp)==1:
        in_dir='./'
        snap = snaphot_file
    else:
        in_dir=''
        for j in range(len(arp)-1): in_dir = in_dir+arp[j]+'/'
        snap = arp[-1]

    print('Generate the files needed by the clustering algorithm...')
    file_midplane = midplane_potential(in_dir, snap, out_dir, align_with=align_with, radius_align=radius_align, halo_id=halo_id, verbose=verbose)
    file_potential = star_potential(in_dir, snap, out_dir, halo_id=halo_id, verbose=verbose)
    GMM_in = GMM_input_file(in_dir, snap, out_dir, file_potential.split('/')[-1], file_midplane.split('/')[-1], halo_id=halo_id, verbose=verbose)

    print('Run the clustering...')
    file_dec = gmm_clustering_for_stars(out_dir, GMM_in.split('/')[-1], number_of_clusters = number_of_clusters,
                             covariance_type = covariance_type, whiten_data = whiten_data, n_init=n_init, plot=plot, verbose=verbose)

    print('Plot the results...')
    draw_2d_distributions(in_dir, snap, out_dir, file_midplane.split('/')[-1], file_dec.split('/')[-1], halo_id=halo_id, verbose=verbose)

    print('Finish')
    return

######################################################################################################################################


def star_potential(in_dir, snap, out_dir, halo_id=1, verbose=True):

    if verbose:
        print('-------------------------------------------------------------------------------------------------------------------------')
        print('This functions calls the Fortran90 parallel module twobody to compute the gravitational potential')
        print('Required input:')
        print('in_dir = the path to the directory where the simulation is')
        print('snap = the filename of the simulation snaphot you want to process')
        print('out_dir = the path to the directory where the output file will be saved')
        print('If you want to run it for other halo than the most massive one in the simulation snapshot, you need to pass it also the halo_id arg')

    start0 = time.time()
    if in_dir[-1]!='/': in_dir = in_dir+'/'
    filein = in_dir+snap
    if out_dir[-1]!='/': out_dir = out_dir+'/'
    fileout = out_dir+snap+('.halo_%i'%halo_id)+'.star_potential.dat'

    if os.path.isfile(fileout):
        print('The potential file already exists. Returning...')
        return fileout

    print('Running twobody.star_potential for snaphot file %s ...'%filein)

    try:
        s = pynbody.load(filein)
    except:
        print('Could not load with pynbody the file %s'%filein)
        print('Check your file name! It should contain the complete path.')
        return

    h = s.halos()
    h1 = h[halo_id]
    s.physical_units()

    cen = pynbody.analysis.halo.center(h1,mode='ssc',retcen=True)
    for j, arr in enumerate(['x','y','z']):
        s[arr] -= cen[j]
    rvir = np.float(pynbody.array.SimArray(np.max(h1['r'].in_units('kpc')),'kpc'))
    mvir = np.sum(h1['mass'].in_units('Msol'))
    vvir = np.sqrt(grav_const*np.float(mvir)/np.float(rvir))

    eps = h1['eps'].min()
    iord = np.array(h1.s['iord'])
    md = np.array(h1.d['mass'].in_units('Msol'))
    xd = np.array(h1.d['x'].in_units('kpc'))
    yd = np.array(h1.d['y'].in_units('kpc'))
    zd = np.array(h1.d['z'].in_units('kpc'))
    nd = len(md)
    mg = np.array(h1.g['mass'].in_units('Msol'))
    xg = np.array(h1.g['x'].in_units('kpc'))
    yg = np.array(h1.g['y'].in_units('kpc'))
    zg = np.array(h1.g['z'].in_units('kpc'))
    ng = len(mg)
    ms = np.array(h1.s['mass'].in_units('Msol'))
    xs = np.array(h1.s['x'].in_units('kpc'))
    ys = np.array(h1.s['y'].in_units('kpc'))
    zs = np.array(h1.s['z'].in_units('kpc'))
    ns = len(ms)
    start1 = time.time()
    phi = twobody.star_potential(md,xd,yd,zd,mg,xg,yg,zg,ms,xs,ys,zs,eps)
    end1 = time.time()

    phi = -phi*grav_const
    f = open(fileout,'w')
    pickle.dump({'snap':snap,'z':s.properties['z'],'time':s.properties['time'].in_units('Gyr'),'halo_id':halo_id,
                'cen':cen,'rvir':rvir,'mvir':mvir,'vvir':vvir,'eps':eps,'nd':nd,'ng':ng,'ns':ns,
                'iord':iord,'phi':phi},f)
    f.close()
    del(f)

    del(iord,ms,xs,ys,zs,mg,xg,yg,zg,md,xd,yd,zd,phi)
    del(h1,h,s)
    end = time.time()
    if verbose:
        print('Time to load the simulation snapshot: %s'%secondsToStr(start1-start0))
        print('Time to compute the potential for stars: %s'%secondsToStr(end1-start1))
        print('Total runtime: %s'%secondsToStr(end-start0))
        print('-------------------------------------------------------------------------------------------------------------------------')
    return fileout

######################################################################################################################################


def midplane_potential(in_dir, snap, out_dir, align_with='baryon', radius_align=None, halo_id=1, verbose=True, *args, **kwargs):

    if verbose:
        print('-----------------------------------------------------------------------------------------------------------------------------------------------------')
        print('This functions calls the Fortran90 parallel module twobody to compute the gravitational potential in the equatorial plane.')
        print('The z-axis perpendicular to this plane is defined by the angular momentum of either stars, gas, or stars+gas inside 0.1rvir. Default uses baryons.')
        print('Set the arg align_with to use the stars or the gas for the alignement.')
        print('Set the arg radius_align to override the 0.1rvir default. The value can be either in fraction of rvir or in kpc.')
        print('Required input:')
        print('in_dir = the path to the directory where the simulation is')
        print('snap = the filename of the simulation snap you want to process')
        print('out_dir = the path to the directory where the output file will be saved')
        print('If you want to run it for other halo than the most massive one in the simulation snapshot, you need to pass it also the halo_id arg')

    if align_with not in ['star','gas','baryon']: align_with='baryon'

    start = time.time()
    if in_dir[-1]!='/': in_dir = in_dir+'/'
    filein = in_dir+snap
    if out_dir[-1]!='/': out_dir = out_dir+'/'
    fileout = out_dir+snap+('.halo_%i'%halo_id)+'.align_with_'+align_with+'.midplane_potential.dat'

    if os.path.isfile(fileout):
        print('The potential file already exists. Returning...')
        return fileout

    s = pynbody.load(filein)
    h = s.halos()
    h1 = h[halo_id]
    s.physical_units()
    rmin = 3.*h1['eps'].min()

    cen = pynbody.analysis.halo.center(h1,mode='ssc',retcen=True)
    for j, arr in enumerate(['x','y','z']):
        s[arr] -= cen[j]
    rvir = np.float(pynbody.array.SimArray(np.max(h1['r'].in_units('kpc')),'kpc'))
    mvir = np.sum(h1['mass'].in_units('Msol'))
    vvir = np.sqrt(grav_const*np.float(mvir)/np.float(rvir))
    vcen=pynbody.analysis.halo.vel_center(h1,retcen=True,cen_size=0.1*rvir)
    for j, arr in enumerate(['vx','vy','vz']):
        s[arr] -= vcen[j]

    if radius_align is None: radius_align = 0.1*rvir
    else:
        if radius_align < 1.0: radius_align = radius_align*rvir

    if align_with=='baryon':
        js1 = pynbody.analysis.angmom.ang_mom_vec(h1[filt.BandPass('r',rmin,radius_align)].g)
        js2 = pynbody.analysis.angmom.ang_mom_vec(h1[filt.BandPass('r',rmin,radius_align)].s)
        js = js1+js2
        print('halo aligned with the angular momentum of the inner baryons')
    if align_with=='star':
        js = pynbody.analysis.angmom.ang_mom_vec(h1[filt.BandPass('r',rmin,radius_align)].s)
        print('halo aligned with the angular momentum of the inner stars')
    if align_with=='gas':
        js = pynbody.analysis.angmom.ang_mom_vec(h1[filt.BandPass('r',rmin,radius_align)].g)
        print('halo aligned with the angular momentum of the inner gas')
    trans = pynbody.analysis.angmom.calc_faceon_matrix(js)
    s.transform(trans)

    eps = h1['eps'].min()
    minR = 0.1*eps
    maxR = rvir
    nbins = 100
    bin_edges = np.logspace(np.log10(minR), np.log10(maxR), num=nbins+1)
    Rbins = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    xin = np.concatenate((Rbins,-Rbins,np.zeros(len(Rbins)),np.zeros(len(Rbins))))
    yin = np.concatenate((np.zeros(len(Rbins)),np.zeros(len(Rbins)),Rbins,-Rbins))
    xin = xin.flatten()
    yin = yin.flatten()
    zin = np.zeros(len(xin))
    ni = len(xin)

    md = np.array(h1.d['mass'].in_units('Msol'))
    xd = np.array(h1.d['x'].in_units('kpc'))
    yd = np.array(h1.d['y'].in_units('kpc'))
    zd = np.array(h1.d['z'].in_units('kpc'))
    nd = len(md)
    mg = np.array(h1.g['mass'].in_units('Msol'))
    xg = np.array(h1.g['x'].in_units('kpc'))
    yg = np.array(h1.g['y'].in_units('kpc'))
    zg = np.array(h1.g['z'].in_units('kpc'))
    ng = len(mg)
    ms = np.array(h1.s['mass'].in_units('Msol'))
    xs = np.array(h1.s['x'].in_units('kpc'))
    ys = np.array(h1.s['y'].in_units('kpc'))
    zs = np.array(h1.s['z'].in_units('kpc'))
    ns = len(ms)

    print('compute the jc-E mapping considering the halo in isolation and recomputing the potential... ')

    start1 = time.time()
    phi = twobody.midplane_potential(md,xd,yd,zd,mg,xg,yg,zg,ms,xs,ys,zs,xin,yin,zin,eps)
    end1 = time.time()

    phi = -phi*grav_const
    pot_midplane = np.zeros(nbins)
    for j in range(nbins):
        pot_midplane[j] = (phi[j]+phi[j+nbins]+phi[j+2*nbins]+phi[j+3*nbins])*0.25

    start2 = time.time()
    vcirc2 = twobody.midplane_vcirc2(md,xd,yd,zd,mg,xg,yg,zg,ms,xs,ys,zs,xin,yin,zin,eps)
    end2 = time.time()

    v_circ = np.zeros(nbins)
    for j in range(nbins):
        v_circ[j] = np.sqrt(grav_const*(vcirc2[j]+vcirc2[j+nbins]+vcirc2[j+2*nbins]+vcirc2[j+3*nbins])*0.25)
    j_circ = v_circ*Rbins
    bindingE = 0.5*(v_circ**2) + pot_midplane

    f = open(fileout,'w')
    pickle.dump({'snap':snap,'z':s.properties['z'],'time':s.properties['time'].in_units('Gyr'),'halo_id':halo_id,
                'cen':cen,'vcen':vcen,'matrix':trans,'rvir':rvir,'mvir':mvir,'vvir':vvir,'radius_align':radius_align,'eps':eps,
                'R':Rbins,'v_circ':v_circ,'j_circ':j_circ,'bindingE':bindingE},f)
    f.close()
    del(f)

    del(h1)
    del(h)
    del(s)
    end = time.time()
    print('Runtime: %s'%secondsToStr(end-start))
    print('-------------------------------------------------------------------------------------------------------------------------')
    return fileout

######################################################################################################################################


def GMM_input_file(in_dir, snap, out_dir, file_potential, file_midplane, halo_id=1, verbose=True):

    if verbose:
        print('-------------------------------------------------------------------------------------------------------------------------')
        print('This function computes the circularity distribution for the stars in one halo of simulation snap.')
        print('By default it will use the most massive halo (halo_id=1).')
        print('The halo is orientated using the rotation matrix defined for the midplane potential.')
        print('Required input:')
        print('in_dir = the path to the directory where the simulation is')
        print('snap = the filename of the simulation snaphot you want to process')
        print('out_dir = the path to the directory where the output file will be saved')
        print('file_potential = the file containing the gravitational potential of all halo stars')
        print('file_midplane = the file containing the definition of the equatorial plane')


    start = time.time()
    if in_dir[-1]!='/': in_dir = in_dir+'/'
    filein = in_dir+snap
    if out_dir[-1]!='/': out_dir = out_dir+'/'

    print('Checking if the halo_id and the name of the file_potential and file_midplane are consistent...')
    arp = file_potential.split('.')
    if ('.halo_%i'%halo_id) not in arp:
        file_potential = file_potential.replace('.halo_*',('.halo_%i'%halo_id),1)
        print ('I will be using as file_potential: %s '%file_potential)
    arp = file_midplane.split('.')
    if ('.halo_%i'%halo_id) not in arp:
        file_midplane = file_midplane.replace('.halo_*',('.halo_%i'%halo_id),1)
        print('I will be using as file_midplane: %s '%file_midplane)
    file_potential = out_dir + file_potential
    file_midplane = out_dir + file_midplane

    fileout = out_dir+snap+('.halo_%i'%halo_id)+'.gmm_input_space.dat'
    print('The output will be saved to %s'%fileout)

    s = pynbody.load(filein)
    h = s.halos()
    h1 = h[halo_id]
    s.physical_units()

    try:
        data = pickle.load(open(file_midplane,'r'))
    except:
        print('Input file %s does not exist! Returning...'%file_midplane)
        return
    trans = data['matrix']
    cen = data['cen']
    vcen = data['vcen']
    rvir = data['rvir']
    mvir = data['mvir']
    vvir = data['vvir']
    j_circ = data['j_circ']
    bindingE = data['bindingE']
    for j, arr in enumerate(['x','y','z']):
        s[arr] -= cen[j]
    for j, arr in enumerate(['vx','vy','vz']):
        s[arr] -= vcen[j]
    s.transform(trans)

    h1['jz'] = h1['x']*h1['vy']-h1['y']*h1['vx']
    h1['jx'] = h1['y']*h1['vz']-h1['z']*h1['vy']
    h1['jy'] = h1['z']*h1['vx']-h1['x']*h1['vz']
    h1['jp'] = np.sqrt(h1['jx']**2+h1['jy']**2)

    mass = (h1.s['mass'].in_units('Msol')).view(np.ndarray)
    jz = (h1.s['jz'].in_units('kpc km s**-1')).view(np.ndarray)
    jp = (h1.s['jp'].in_units('kpc km s**-1')).view(np.ndarray)
    ke = (h1.s['ke'].in_units('km**2 s**-2')).view(np.ndarray)
    iord = np.array(h1.s['iord'])
    srt = np.argsort(iord)
    mass = mass[srt]
    jz = jz[srt]
    jp = jp[srt]
    ke = ke[srt]
    iord = iord[srt]

    try:
        data = pickle.load(open(file_potential,'r'))
    except:
        print('Input file %s does not exist! Returning...'%file_potential)
        return
    iord_pot = data['iord']
    phi = data['phi']
    srt = np.argsort(iord_pot)
    iord_pot = iord_pot[srt]
    phi = phi[srt]

    print('check if the two iord arrays are the same...')
    if len(iord)!=len(iord_pot):
        print('the two iord arrays should have the same length! returning')
        return

    if abs(np.sum(iord-iord_pot))>0:
        print('the two iord arrays are not equal! returning')
        return

    energy = ke + phi

    reduce_e = 10.**(-int(np.log10(max(abs(bindingE)))))
    reduce_jc = 10.**(-int(np.log10(max(j_circ))))
    interp = scipy.interpolate.InterpolatedUnivariateSpline(bindingE*reduce_e, j_circ*reduce_jc, k=1)
    rep = interp(energy*reduce_e)
    jc = rep/reduce_jc

    f = open(fileout,'w')
    pickle.dump({'snap':snap,'z':s.properties['z'],'time':s.properties['time'].in_units('Gyr'),'halo_id':halo_id,
                'cen':cen,'vcen':vcen,'matrix':trans,'rvir':rvir,'mvir':mvir,'vvir':vvir,
                'mass':mass,'jz':jz,'jp':jp,'jc':jc,'energy':energy,'iord':iord},f)
    f.close()
    del(f)

    return fileout

######################################################################################################################################


def GMM_input_space(out_dir, GMM_input_file, verbose=True):

    if verbose:
        print('-------------------------------------------------------------------------------------------------------------------------')
        print('This function prepares the input array for the clustering computes.')
        print('Required input:')
        print('out_dir = the path to the directory where the GMM_input_file is')
        print('GMM_input_file = the file containing the normalized angular momenta and binding energies of stars')


    if out_dir[-1]!='/': out_dir = out_dir+'/'
    GMM_input_file = out_dir+GMM_input_file

    try:
        data = pickle.load(open(GMM_input_file,'r'))
    except:
        print('Input file %s  does not exist! Returning...'%GMM_input_file)
        return np.array(1)*np.nan, np.array(1)*np.nan

    iord = data['iord']
    mass = data['mass']
    ns = len(iord)
    obs = np.ndarray(shape=(ns,3), dtype=float)
    obs[:,0] = data['jz']/data['jc']
    obs[:,1] = data['jp']/data['jc']
    obs[:,2] = data['energy']/max(abs(data['energy']))

    initial_length = ns

    varr = obs[:,0].flatten()
    print('The values of jz/jc should be within [-1,1], actual limits: %.2f, %.2f'%(min(varr),max(varr)))
    keep = np.where((varr >= -1.5) & (varr <= 1.5))
    print('I will keep only entries within [-1.5,1.5], so i need to discard this fraction: %.4f'%(1.-float(len(varr[keep])/float(initial_length))))
    obs = obs[keep]
    iord = iord[keep]
    mass = mass[keep]

    varr = obs[:,1].flatten()
    print('The values of jp/jc should be within [0,1], actual limits: %.2f, %.2f'%(min(varr),max(varr)))
    keep = np.where((varr >= 0) & (varr <= 1.5))
    print('I will keep only entries within [0,1.5], so i need to discard this fraction: %.4f'%(1.-float(len(varr[keep])/float(initial_length))))
    obs = obs[keep]
    iord = iord[keep]
    mass = mass[keep]

    varr = obs[:,2].flatten()
    print('The values of normalized binding energy, e, should be within [-1,0], actual limits: %.2f, %.2f'%(min(varr),max(varr)))
    keep = np.where(varr <= 0.)
    print('I will keep only entries within [-1,0], so i need to discard this fraction: %.4f'%(1.-float(len(varr[keep])/float(initial_length))))
    obs = obs[keep]
    iord = iord[keep]
    mass = mass[keep]

    return obs, iord, mass

###############################################################################################


def gmm_clustering_for_stars(out_dir, GMM_input_file, n_init=1, plot=True, verbose=True, *args, **kwargs):

    if verbose:
        print('-------------------------------------------------------------------------------------------------------------------------')
        print('This function calls scikit learn to run gaussian mixture models.')
        print('Required input:')
        print('out_dir = the path to the directory, where the GMM_input_file is')
        print('GMM_input_file = the file containing the normalized angular momenta and binding energies of stars')
        print('By default it will require the algorithm to search for 2 clusters.')
        print('If you want to run it with a different number of clusters, set the arg number_of_clusters.')
        print('By default it assumes each clustering component has it own general covariance matrix, arg covariance_type=full')
        print('By default the input data is not centered to the mean and scaled to unit variance. Set arg whiten_data=True to override the default.')
        print('By default, the Expectation-Maximization algorithm is run for up to 100 times or up until the convergence threshold is met (1.0e-03)')
        print('using a single k-means initialization (n_init=1).')
        print('By default it will also plot the mass weight 1D distribution of the clusters found in the original input space of jz/jc, jp/jc and binding energy.')
        print('Set arg plot=False if these plots are not to be saved.')
        print('If you want to use other args for the clustering algorithm, check the Gaussian Mixture webpage at http://scikit-learn.org/stable/')

    number_of_clusters = kwargs.get('number_of_clusters', None)
    whiten_data = kwargs.get('whiten_data', None)
    covariance_type = kwargs.get('covariance_type', None)

    obs, iord, mass = GMM_input_space(out_dir,GMM_input_file,verbose=verbose)
    if np.isnan(iord[0]):
        print('Could not get the input space matrix. Returning...')
        return

    nf = 3
    ns = len(obs[:,0].flatten())

    if out_dir[-1]!='/': out_dir = out_dir+'/'
    filename_base = out_dir+GMM_input_file.replace('gmm_input_space.dat','inputspace_jzjc_jpjc_energy',1)
    name_add = '_'
    white = ''

    if whiten_data is not None:
        inputdata = scale(obs)
        print('Centering data to the mean and component wise scale to unit variance. Output file names contain _white')
        white = '_white'
    else: inputdata = obs

    if number_of_clusters is None: number_of_clusters=2
    print('searching for %i clusters'%number_of_clusters)


    print('******************************************************************************************')
    print('Running Gaussian Mixture Model clustering in the space of jz/jc, jp/jc, binding energy')

    if covariance_type not in ['full','tied','diag','spherical']: covariance_type='full'

    filename_out=filename_base+'.scikit_gmm_'+covariance_type+'_'+("%i" % number_of_clusters)+'clusters'+white

    st = time.time()
    aclus = GMM(n_components=number_of_clusters, covariance_type=covariance_type, n_init=n_init)
    aclus.fit(inputdata)
    centre = aclus.means_ # The means of the adjusted gaussians
    weight = aclus.weights_ # The mixing weights for each mixture component
    covar = aclus.covars_ # Covariance parameters for each mixture component. The shape depends on covariance_type.
    label = aclus.predict(inputdata) # predicted labels
    logprob, posterior_probability = aclus.score_samples(inputdata)
    bic = aclus.bic(inputdata)
    elapsed_time = (time.time() - st)/60.
    print("Elapsed time running GMM: %.2f min" % elapsed_time)

    f = open(filename_out+'.dat','w')
    pickle.dump({'input_file':GMM_input_file,'runtime_min':elapsed_time,
                'gmeans':centre,'gweight':weight,'gcovar':covar,'bic':bic,
                'label':label,'p_label':posterior_probability,'iord':iord},f)
    f.close()
    del(f)

    print('number of clusters: %i'%number_of_clusters)
    print('means:              ',centre)
    print('weights:            ',weight)
    print('covariance matrix:  ',covar)
    print('bic:                ',bic)
    print('Save GMM clustering result to %s.dat'%filename_out)

    if plot:
        print('Plotting the results as mass weighted distributions in the original feature space of jz/jc, jp/jc, energy')
        indices = np.unique(label)
        indices = indices[np.argsort(indices)]
        nk = len(indices)
        jzjc = np.array(obs[:,0]).flatten()
        jpjc = np.array(obs[:,1]).flatten()
        energy = np.array(obs[:,2]).flatten()
        nbin_jzjc = 301
        range_jzjc = (-1.5,1.5)
        nbin_jpjc = 150
        range_jpjc = (0,1.5)
        nbin_energy = 100
        range_energy = (-1.0,1)

        mass_total_jzjc, bin_edges, binnumber = scipy.stats.binned_statistic(jzjc, mass, statistic='sum', bins=nbin_jzjc, range=range_jzjc)
        jzjc_bins = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        mass_total_jpjc, bin_edges, binnumber = scipy.stats.binned_statistic(jpjc, mass, statistic='sum', bins=nbin_jpjc, range=range_jpjc)
        jpjc_bins = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        mass_total_energy, bin_edges, binnumber = scipy.stats.binned_statistic(energy, mass, statistic='sum', bins=nbin_energy, range=range_energy)
        energy_bins = 0.5 * (bin_edges[:-1] + bin_edges[1:])

        ymax = np.nanmax(np.array([np.nanmax(mass_total_jzjc),np.nanmax(mass_total_jpjc),np.nanmax(mass_total_energy)]))

        color_scheme = plt.get_cmap('jet_r')
        cNorm  = colors.Normalize(vmin=0, vmax=1)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=color_scheme)
        colorzs = []
        deltac = 1./(nk-1)
        for k in range(nk): colorzs.append(scalarMap.to_rgba(k*deltac))

        plt.close()
        fig, axs = plt.subplots(1,3, figsize=(10.5,3.5))
        fig.subplots_adjust(left=0.05,bottom=0.15,right=0.97,top=0.93,hspace=0.13,wspace=0.13)
        for ax, j in zip(axs.ravel(),[0,1,2]):
            plt.setp(ax.get_xticklabels(),fontsize=12)
            plt.setp(ax.get_yticklabels(),fontsize=12)
            ax.xaxis.labelpad = 1
            ax.yaxis.labelpad = 1
            #ax.set_ylim(0,1.1*ymax)
            #if j in [1,2]: plt.setp(ax.get_yticklabels(), visible=False)
            if j==0:
                ax.set_xlim(-1.5,1.5)
                ax.set_xlabel(r"j$_{\rm z}$/j$_{\rm c}$",fontsize=18)
                ax.set_ylabel(r"M$_{\rm *}$ [M$_{\rm\odot}$]")
                ax.plot(jzjc_bins,mass_total_jzjc,ls='-',color='grey',lw=2)
                x = jzjc
                rangep = range_jzjc
                xbins = jzjc_bins
                bins = nbin_jzjc
                text_x_pos = 0.03
                text_x_align = 'left'
            if j==1:
                ax.set_xlim(0,1.5)
                ax.set_xlabel(r"j$_{\rm p}$/j$_{\rm c}$",fontsize=18)
                ax.plot(jpjc_bins,mass_total_jpjc,ls='-',color='grey',lw=2)
                x = jpjc
                rangep = range_jpjc
                xbins = jpjc_bins
                bins = nbin_jpjc
                text_x_pos = 0.97
                text_x_align = 'right'
            if j==2:
                ax.set_xlim(-1,0)
                ax.set_xlabel(r"e/max(|e|)",fontsize=18)
                ax.plot(energy_bins,mass_total_energy,ls='-',color='grey',lw=2)
                x = energy
                rangep = range_energy
                xbins = energy_bins
                bins = nbin_energy
                text_x_pos = 0.97
                text_x_align = 'right'
            for k,indx in zip(range(nk),indices):
                massink, bin_edges, binnumber = scipy.stats.binned_statistic(x[label==indx], mass[label==indx], statistic='sum', bins=bins, range=rangep)
                massk = np.multiply(mass,np.ravel(np.array(posterior_probability[:,indx]).flatten()))
                massinkp, bin_edges, binnumber = scipy.stats.binned_statistic(x, massk, statistic='sum', bins=bins, range=rangep)
                ax.plot(xbins, massink, color=colorzs[k], ls='-', lw=1.5)
                ax.plot(xbins, massinkp, color=colorzs[k], ls='--', lw=1.5)
                strt = r"$\rm<x_{\rm %i}>$=%.2f"%(indx,np.sum(x*massk)/np.sum(massk))
                ax.text(text_x_pos,0.90-k*0.08,strt,color=colorzs[k],fontsize=15,horizontalalignment=text_x_align,verticalalignment='center',transform = ax.transAxes)
        plt.savefig(filename_out+'.png')
        plt.close()

    return (filename_out+'.dat')

###############################################################################################


def ellipticity_from_moments(x,y,weights):
    r = np.sqrt(x**2+y**2)
    f = weights
    x = x[np.isfinite(f)]
    y = y[np.isfinite(f)]
    r = r[np.isfinite(f)]
    f = f[np.isfinite(f)]
    x = x[r>0]
    y = y[r>0]
    f = f[r>0]
    r = r[r>0]
    Mxx = np.sum(f*(x/r)**2)/np.sum(f)
    Myy = np.sum(f*(y/r)**2)/np.sum(f)
    Mxy = np.sum(f*(x*y)/r**2)/np.sum(f)
    Q = Mxx - Myy # = a-b/a+b cos(2phi)
    U = 2.*Mxy    # = a-b/a+b sin(2phi)
    assymetry = (1.-np.sqrt(Q**2+U**2))/(1.+np.sqrt(Q**2+U**2))
    ellipticity = 1.-assymetry
    return ellipticity

###############################################################################################


def draw_2d_distributions(in_dir, snap, out_dir, file_midplane, file_dec, halo_id=1, verbose=True, *args, **kwargs):

    if verbose:
        print('-------------------------------------------------------------------------------------------------------------------------')
        print('This function will create for each component found in file_dec a png figure with ')
        print('the face-on (left) and edge-on (center) surface mass density, and the edge-on mass weighted line-of-sight velocity (right)')
        print('Required input:')
        print('in_dir = the path to the directory where the simulation is')
        print('snap = the filename of the simulation snap you want to process')
        print('out_dir = the path to the directory where the output file will be saved')
        print('file_midplane = the file containing the definition of the equatorial plane')
        print('file_dec = the file containing the result of the clustering algorithm for the stars of halo halo_id (default halo_id=1)')
        print('The two optional arguments are the number of pixels per side nxny, and the dimension of the field of view fov. ')

    nxny = kwargs.get('nxny', None) # number of pixels per side
    fov = kwargs.get('fov', None) # figure size in kpc. Defaults to 0.1rvir

    start = time.time()
    if in_dir[-1]!='/': in_dir = in_dir+'/'
    filein = in_dir+snap
    if out_dir[-1]!='/': out_dir = out_dir+'/'

    print('Checking if the halo_id and the name of the file_midplane are consistent...')
    arp = file_midplane.split('.')
    if ('.halo_%i'%halo_id) not in arp:
        file_midplane = file_midplane.replace('.halo_*',('.halo_%i'%halo_id),1)
        print ('I will be using as file_midplane: %s '%file_midplane)
    file_midplane = out_dir + file_midplane
    try:
        data_midplane = pickle.load(open(file_midplane,'r'))
    except:
        print('Input file  does not exist! Returning...'%file_midplane)
        return
    cen = data_midplane['cen']
    vcen = data_midplane['vcen']
    trans = data_midplane['matrix']
    rvir = data_midplane['rvir']

    s = pynbody.load(filein)
    h = s.halos()
    h1 = h[halo_id]
    s.physical_units()
    for j, arr in enumerate(['x','y','z']): s[arr] -= cen[j]
    for j, arr in enumerate(['vx','vy','vz']): s[arr] -= vcen[j]
    s.transform(trans)

    file_dec = out_dir + file_dec
    try:
        data_dec = pickle.load(open(file_dec,'r'))
    except:
        print('Input file does not exist. You need to run gmm_clustering_for_stars first!')
        return
    figure_name = file_dec[:-3]+'xy_and_yz_projections.png'

    label = data_dec['label']
    plabel = data_dec['p_label']
    iord_dec = data_dec['iord']
    srt = np.argsort(iord_dec)
    label = label[srt]
    plabel = plabel[srt]
    iord_dec = iord_dec[srt]

    h1.s = h1.s[np.in1d(h1.s['iord'],iord_dec)]
    iord = np.array(h1.s['iord'])
    srt = np.argsort(iord)
    mass = np.array(h1.s['mass'].in_units('Msol'))[srt]
    x = np.array(h1.s['x'].in_units('kpc'))[srt]
    y = np.array(h1.s['y'].in_units('kpc'))[srt]
    z = np.array(h1.s['z'].in_units('kpc'))[srt]
    vx = np.array(h1.s['vx'])[srt]
    iord = iord[srt]

    if fov is None: fov = 0.2*rvir
    lim = 0.5*fov
    if nxny is None:
        pixel_size_kpc = 0.3
        nxny = int(2*fov/pixel_size_kpc)
    rangex = (-lim,lim)
    rangey = (-lim,lim)
    range2d = (rangex,rangey)
    extent=[rangex[0],rangex[1],rangey[0],rangey[1]]

    color_scheme = plt.get_cmap('jet_r')
    cNorm  = colors.Normalize(vmin=0, vmax=1)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=color_scheme)

    total_mass = np.sum(mass)
    indices = np.unique(label)
    indices = indices[np.argsort(indices)]
    nk = len(indices)
    deltac = 1./(nk-1)

    kell = np.zeros(nk)
    kname = []
    for k,indx in zip(range(nk),indices):
        ell = ellipticity_from_moments(y[label==indx],z[label==indx],mass[label==indx])
        kell[k] = ell
        kname.append('GMM component %i'%indx)
    if nk==2:
        if kell[0]<kell[1]:
            kname[0] = 'GMM spheroid'
            kname[1] = 'GMM disk'
            print('The cluster 0 is the spheroid and 1 is the disk')
        else:
            kname[1] = 'GMM spheroid'
            kname[0] = 'GMM disk'
            print('The cluster 1 is the spheroid and 0 is the disk')

    for k,indx in zip(range(nk),indices):
        c_x = x[label==indx]
        c_y = y[label==indx]
        c_z = z[label==indx]
        c_v = vx[label==indx]
        c_mass = mass[label==indx]
        color = scalarMap.to_rgba(k*deltac)
        faceon_mass, xedges, yedges, binnumber = scipy.stats.binned_statistic_2d(c_x, c_y, c_mass, statistic='sum', bins=nxny, range=range2d)
        edgeon_mass, xedges, yedges, binnumber = scipy.stats.binned_statistic_2d(c_y, c_z, c_mass, statistic='sum', bins=nxny, range=range2d)
        edgeon_vlos, xedges, yedges, binnumber = scipy.stats.binned_statistic_2d(c_y, c_z, c_mass*c_v, statistic='sum', bins=nxny, range=range2d)
        edgeon_vlos = edgeon_vlos/edgeon_mass # units = km/s
        faceon_surfacemassdensity = faceon_mass/(pixel_size_kpc*1.0e6) # units = Msol/pc^2
        edgeon_surfacemassdensity = edgeon_mass/(pixel_size_kpc*1.0e6) # units = Msol/pc^2
        mass_fraction = np.sum(c_mass)/total_mass
        plot_name = figure_name[:-3]+('component_%i.png'%indx)
        plot_logmu_vlos(faceon_surfacemassdensity,edgeon_surfacemassdensity,edgeon_vlos,mass_fraction,kname[k],extent,rangex,rangey,color,plot_name)
    color = 'grey'
    faceon_mass, xedges, yedges, binnumber = scipy.stats.binned_statistic_2d(x, y, mass, statistic='sum', bins=nxny, range=range2d)
    edgeon_mass, xedges, yedges, binnumber = scipy.stats.binned_statistic_2d(y, z, mass, statistic='sum', bins=nxny, range=range2d)
    edgeon_vlos, xedges, yedges, binnumber = scipy.stats.binned_statistic_2d(y, z, mass*vx, statistic='sum', bins=nxny, range=range2d)
    edgeon_vlos = edgeon_vlos/edgeon_mass # units = km/s
    faceon_surfacemassdensity = faceon_mass/(pixel_size_kpc*1.0e6) # units = Msol/pc^2
    edgeon_surfacemassdensity = edgeon_mass/(pixel_size_kpc*1.0e6) # units = Msol/pc^2
    mass_fraction = 1.
    plot_name = figure_name[:-3]+'component_all.png'
    plot_logmu_vlos(faceon_surfacemassdensity,edgeon_surfacemassdensity,edgeon_vlos,mass_fraction,'All stars',extent,rangex,rangey,color,plot_name)

    return

###############################################################################################


def plot_logmu_vlos(faceon_surfacemassdensity,edgeon_surfacemassdensity,edgeon_vlos,mass_fraction,figure_title,extent,rangex,rangey,color,plot_name):

    smd_color = "Blues"
    vlos_color = plt.get_cmap('jet')

    plt.close()
    fig = plt.figure(figsize=(10.6,3.5))
    gs = gridspec.GridSpec(1, 3)

    ax1 = plt.subplot(gs[0])
    plt.setp(ax1.get_xticklabels(), fontsize=10)
    plt.setp(ax1.get_yticklabels(), fontsize=10)
    ax1.set_xlabel(r"[kpc]", fontsize=10)
    ax1.set_ylabel(r"[kpc]", fontsize=10)
    ax1.xaxis.labelpad = 2
    ax1.yaxis.labelpad = 2
    ax1.set_xlim(rangex)
    ax1.set_ylim(rangey)
    vmin = np.log10(np.min(faceon_surfacemassdensity[faceon_surfacemassdensity>0.]))
    vmax = np.log10(np.max(faceon_surfacemassdensity[faceon_surfacemassdensity>0.]))
    arr = np.log10(faceon_surfacemassdensity)
    im1 = ax1.imshow(arr.transpose(), interpolation='nearest',extent=extent, aspect='auto', origin='lower',cmap=smd_color, vmin=vmin, vmax=vmax, zorder=1)
    cax1 = fig.add_axes([0.306, 0.15, 0.01, 0.78])
    cb1 = fig.colorbar(im1, cax=cax1,orientation="vertical")
    cb1.set_label(label=r"log($\rm\Sigma/M_{\rm\odot}pc^{\rm -2}$)",fontsize=10,labelpad=-5)
    cb1.ax.tick_params(axis='y',labelleft='off',direction='in',labelright='on',pad=1,size=2)
    for t in cb1.ax.get_yticklabels(): t.set_fontsize(8)

    ax2 = plt.subplot(gs[1])
    plt.setp(ax2.get_xticklabels(), fontsize=10)
    plt.setp(ax2.get_yticklabels(), visible=False)
    ax2.set_xlabel(r"[kpc]", fontsize=10)
    ax2.xaxis.labelpad = 2
    ax2.yaxis.labelpad = 2
    ax2.set_xlim(rangex)
    ax2.set_ylim(rangey)
    vmin = np.log10(np.min(edgeon_surfacemassdensity[edgeon_surfacemassdensity>0.]))
    vmax = np.log10(np.max(edgeon_surfacemassdensity[edgeon_surfacemassdensity>0.]))
    arr = np.log10(edgeon_surfacemassdensity)
    sb = arr
    sb[np.isnan(arr)] = -100.
    im2 = ax2.imshow(arr.transpose(), interpolation='nearest',extent=extent, aspect='auto', origin='lower',cmap=smd_color, vmin=vmin, vmax=vmax, zorder=1)
    cax2 = fig.add_axes([0.618, 0.15, 0.01, 0.78])
    cb2 = fig.colorbar(im2, cax=cax2,orientation="vertical")
    cb2.set_label(label=r"log($\rm\Sigma/M_{\rm\odot}pc^{\rm -2}$)",fontsize=10,labelpad=-5)
    cb2.ax.tick_params(axis='y',labelleft='off',direction='in',labelright='on',pad=1,size=2)
    for t in cb2.ax.get_yticklabels(): t.set_fontsize(8)

    ax3 = plt.subplot(gs[2])
    plt.setp(ax3.get_xticklabels(), fontsize=10)
    plt.setp(ax3.get_yticklabels(), visible=False)
    ax3.set_xlabel(r"[kpc]", fontsize=10)
    ax3.xaxis.labelpad = 2
    ax3.yaxis.labelpad = 2
    ax3.set_xlim(rangex)
    ax3.set_ylim(rangey)
    arr = edgeon_vlos
    sb_min = np.min(sb[sb>-99])
    sb_max = np.max(sb[sb>-99])
    sb_limit = sb_min+0.3*(sb_max-sb_min)
    arr[sb < sb_limit] = np.nan
    vmax = np.nanmax(abs(arr))
    vmin = -vmax
    im3 = ax3.imshow(arr.transpose(), interpolation='nearest',extent=extent, aspect='auto', origin='lower',cmap=vlos_color, vmin=vmin, vmax=vmax, zorder=1)
    cax3 = fig.add_axes([0.93, 0.15, 0.01, 0.78])
    cb3 = fig.colorbar(im3, cax=cax3,orientation="vertical")
    cb3.set_label(label=r"v$_{\rm los}$ [km/s]",fontsize=10,labelpad=-5)
    cb3.ax.tick_params(axis='y',labelleft='off',direction='in',labelright='on',pad=1,size=2)
    for t in cb3.ax.get_yticklabels(): t.set_fontsize(8)

    ax2.text(0.5,1.01,figure_title,color=color,fontsize=12,horizontalalignment='center',verticalalignment='bottom',transform=ax2.transAxes)
    gs.update(left=0.05,bottom=0.15,right=0.93, top=0.93, hspace=0.22, wspace=0.22)
    plt.savefig(plot_name)
    plt.close()

    return

###############################################################################################
