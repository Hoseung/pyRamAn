#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>

//using namespace std;

struct Meta{
    int nbodies, halnum, subnum;
    float massp, aexp, omegat, age;
};

// First try, array of a structure.
struct Data{
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
};

void fortran_skip(std::fstream&);

void fortran_read_int(std::fstream&, int &, int);

void fortran_read_float(std::fstream&, float &, int);

void load_meta(std::fstream&,  Meta&);

void load_data(std::fstream&,  Data&, int);

void allocate_data(Data&, int );

void load(std::string&, Meta& haloinfo, Data& halodata);

