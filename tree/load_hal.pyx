# cdef function is local to this source,
# but def function is callable from python.
# So rd_h is the interface between python and functions from the  c++ code c_rd_halo.h.
cimport cython
from libcpp.string cimport string
import numpy as np

from cython.view cimport array as cvarray

cdef extern from "c_rd_halo.h":
    void load(string, Meta&, Data&, int, int)
    void load_d(string, Meta2&, Data2&, int, int)
    cppclass Meta:
        int nbodies, halnum, subnum;
        float massp, aexp, omegat, age;
        int * allID;
    cppclass Meta2:
        int nbodies, halnum, subnum;
        double massp, aexp, omegat, age;
        int * allID;
    cppclass Data:
        int * np;
        int * hnu;
        int * hhost ;
        float * ang ;
        float * mass;
        float * sp;
        float * pos;
        float * vel;
        float * radius;
        float * energy;
        float * vir;
        float * profile;
        float * gal;
        int * g_nbin;
        float * g_rr;
        float * g_rho;
    cppclass Data2:
        int * np;
        int * hnu;
        int * hhost ;
        double * ang ;
        double * mass;
        double * sp;
        double * pos;
        double * vel;
        double * radius;
        double * energy;
        double * vir;
        double * profile;
        double * gal;
        int * g_nbin;
        double * g_rr;
        double * g_rho;

cpdef read_file(filename, nbodies, is_gal):
    cdef Meta haloinfo
    cdef Data halodata
    if is_gal > 0:
        load(filename, haloinfo, halodata, nbodies, 1)
    else:
        load(filename, haloinfo, halodata, nbodies, 0)

    # I don't need to do this.
    # make a wrapper class and pass it directly to python
    # or, pass the pointer to the array to numpy directly (without copying)
    # but.. I don't know how to do that :(
    ntot = haloinfo.halnum + haloinfo.subnum
    allid = np.asarray(<int[:nbodies]> haloinfo.allID)
    hnu = np.asarray(<int[:ntot]> halodata.hnu)
    nump = np.asarray(<int[:ntot]> halodata.np)
    hhost = np.asarray(<int[:ntot*5]> halodata.hhost)
    ang = np.asarray(<float[:ntot*3]> halodata.ang)
    energy = np.asarray(<float[:ntot*3]> halodata.energy)
    mass = np.asarray(<float[:ntot]> halodata.mass)
    vir = np.asarray(<float[:ntot*4]> halodata.vir)
    profile = np.asarray(<float[:ntot*2]> halodata.profile)
    sp = np.asarray(<float[:ntot]> halodata.sp)
    pos = np.asarray(<float[:ntot*3]> halodata.pos)
    vel = np.asarray(<float[:ntot*3]> halodata.vel)
    radius = np.asarray(<float[:ntot*4]> halodata.radius)
    if is_gal > 0:
        gal = np.asarray(<float[:ntot*3]> halodata.gal)
        g_nbin = np.asarray(<int[:ntot]> halodata.g_nbin)
        g_rr = np.asarray(<float[:ntot*100]> halodata.g_rr)
        g_rho = np.asarray(<float[:ntot*100]> halodata.g_rho)
        return [allid, haloinfo.nbodies, haloinfo.halnum, haloinfo.subnum, haloinfo.massp,
            haloinfo.aexp, haloinfo.omegat, haloinfo.age,
            nump, hnu, hhost, ang, energy, mass, radius, pos, sp, vel, vir, profile,
            gal, g_nbin, g_rr, g_rho]
    else:
        return [allid, haloinfo.nbodies, haloinfo.halnum, haloinfo.subnum, haloinfo.massp,
            haloinfo.aexp, haloinfo.omegat, haloinfo.age,
            nump, hnu, hhost, ang, energy, mass, radius, pos, sp, vel, vir, profile]


cpdef read_file_double(filename, nbodies, is_gal):
    cdef Meta2 haloinfo
    cdef Data2 halodata
    if is_gal > 0:
        load_d(filename, haloinfo, halodata, nbodies, 1)
    else:
        load_d(filename, haloinfo, halodata, nbodies, 0)

    # I don't need to do this.
    # make a wrapper class and pass it directly to python
    # or, pass the pointer to the array to numpy directly (without copying)
    # but.. I don't know how to do that :(
    ntot = haloinfo.halnum + haloinfo.subnum
    allid = np.asarray(<int[:nbodies]> haloinfo.allID)
    hnu = np.asarray(<int[:ntot]> halodata.hnu)
    nump = np.asarray(<int[:ntot]> halodata.np)
    hhost = np.asarray(<int[:ntot*5]> halodata.hhost)
    ang = np.asarray(<double[:ntot*3]> halodata.ang)
    energy = np.asarray(<double[:ntot*3]> halodata.energy)
    mass = np.asarray(<double[:ntot]> halodata.mass)
    vir = np.asarray(<double[:ntot*4]> halodata.vir)
    profile = np.asarray(<double[:ntot*2]> halodata.profile)
    sp = np.asarray(<double[:ntot]> halodata.sp)
    pos = np.asarray(<double[:ntot*3]> halodata.pos)
    vel = np.asarray(<double[:ntot*3]> halodata.vel)
    radius = np.asarray(<double[:ntot*4]> halodata.radius)
    if is_gal > 0:
        gal = np.asarray(<double[:ntot*3]> halodata.gal)
        g_nbin = np.asarray(<int[:ntot]> halodata.g_nbin)
        g_rr = np.asarray(<double[:ntot*100]> halodata.g_rr)
        g_rho = np.asarray(<double[:ntot*100]> halodata.g_rho)
        return [allid, haloinfo.nbodies, haloinfo.halnum, haloinfo.subnum, haloinfo.massp,
            haloinfo.aexp, haloinfo.omegat, haloinfo.age,
            nump, hnu, hhost, ang, energy, mass, radius, pos, sp, vel, vir, profile,
            gal, g_nbin, g_rr, g_rho]
    else:
        return [allid, haloinfo.nbodies, haloinfo.halnum, haloinfo.subnum, haloinfo.massp,
            haloinfo.aexp, haloinfo.omegat, haloinfo.age,
            nump, hnu, hhost, ang, energy, mass, radius, pos, sp, vel, vir, profile]
