from libcpp.string cimport string

cdef extern from "c_rd_halo.h":
    int voidfunc()
    void load(string)
    cdef cppclass Meta:
        int nbodies, halnum, subnum;
        float massp, aexp, omegat, age;
    cdef cppclass Data2:
        int * np;
        int * hnu;
        int * hhost ;
        float * ang ;
        float * mass;
        float * mvir;
        float * rvir;
        float * sp;
        float * pos;
        float * vel;
        float * radius;
        float * energy;

def rd_h(fname):
    cdef Meta* haloinfo
    cdef Data2* halodata
    a = voidfunc()
    load(fname)

