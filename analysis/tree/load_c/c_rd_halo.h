#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>

//using namespace std;

struct Meta{
    int nbodies, halnum, subnum;
    float massp, aexp, omegat, age;
};

struct Data2{
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
};

void fortran_skip(std::fstream&);

void fortran_read_int(std::fstream&, int &, int);

void fortran_read_float(std::fstream&, float &, int);

void load_meta(std::fstream&,  Meta&);

void load_data(std::fstream&,  Data2&, int);
void load_data_gal(std::fstream&,  Data2&, int, int);

void allocate_data(Data2&, int, int);

void load(std::string&, Meta&, Data2&, int);
