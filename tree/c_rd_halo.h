#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>

//using namespace std;
struct Meta{
    int nbodies, halnum, subnum;
    float massp, aexp, omegat, age;
    int * allID;
};

struct Meta2{
    int nbodies, halnum, subnum;
    double massp, aexp, omegat, age;
    int * allID;
};

struct Data{
    int * np;
    int * hnu;
    int * hhost;
    long * imbp;
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

struct Data2{
    int * np;
    int * hnu;
    int * hhost;
    long * imbp;
    double * ang;
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
};

/*
Using class may allow me to choose the data type of the member.
*/

void fortran_skip(std::fstream&);
void fortran_read_int(std::fstream&, int &, int);
void fortran_read_long(std::fstream&, long &, int);
void fortran_read_float(std::fstream&, float &, int);
void fortran_read_double(std::fstream&, double &, int);

void load_meta(std::fstream&, Meta&);
void load_meta_d(std::fstream&, Meta2&);

void load_data_gal(std::fstream&,  Data&, int, int, int);
void load_data_gal_d(std::fstream&,  Data2&, int, int, int);

void allocate_data(Data&, int, int, int);
void allocate_data_d(Data2&, int, int, int);

void load(std::string&, Meta&, Data&, int, int, int);
void load_d(std::string&, Meta2&, Data2&, int, int, int);
