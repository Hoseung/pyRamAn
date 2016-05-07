#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include "c_rd_halo.h"

//using namespace std;
// First try, array of a structure.
void fortran_skip(std::fstream& fin){
    int dummy;
    fin.read((char *) &dummy, sizeof(dummy));
    fin.ignore(dummy);
    fin.read((char *) &dummy, sizeof(dummy));
    return ;
}

void fortran_read_int(std::fstream& fin, int &data, int len){
    int dummy;
    fin.read((char *) &dummy, sizeof(dummy));
    fin.read((char *) &data, sizeof(int) * len);
//    std::cout << data << std::endl;
    fin.read((char *) &dummy, sizeof(dummy));
    return ;
}

void fortran_read_float(std::fstream& fin, float &data, int len){
    int dummy;
//    double data;
    fin.read((char *) &dummy, sizeof(dummy));
    fin.read((char *) &data, sizeof(float) * len);
//    std::cout << data << std::endl;
    fin.read((char *) &dummy, sizeof(dummy));
    return ;
}

void load_meta(std::fstream& fin,  Meta& haloinfo){
    int dummy [2];
    fortran_read_int(fin, haloinfo.nbodies, 1);
    fortran_read_float(fin, haloinfo.massp, 1);
    fortran_read_float(fin, haloinfo.aexp, 1);
    fortran_read_float(fin, haloinfo.omegat, 1);
    fortran_read_float(fin, haloinfo.age, 1);
    fortran_read_int(fin, dummy[0], 2);
    haloinfo.halnum = dummy[0];
    haloinfo.subnum = dummy[1];
//    std::cout << " Meta done" << std::endl;
}


void load_data(std::fstream& fin,  Data& halodata, int tothal){
    int* ids;
    for (int i=0; i < tothal; i++){
        fortran_read_int(fin, halodata.np[i], 1);
        std::cout << halodata.np[i] << std::endl;
        ids = new int [halodata.np[i]];
        fortran_read_int(fin, ids[0], halodata.np[i]);
        //delete [] ids;
        fortran_read_int(fin, halodata.hnu[i], 1);
        fortran_skip(fin);// time step
        fortran_read_int(fin, halodata.hhost[i*5], 5);
        fortran_read_float(fin, halodata.mass[i], 1);
        fortran_read_float(fin, halodata.pos[i*3], 3);
        fortran_read_float(fin, halodata.vel[i*3], 3);
        fortran_read_float(fin, halodata.ang[i*3], 3);
        fortran_read_float(fin, halodata.radius[i*4], 4);
        fortran_read_float(fin, halodata.energy[i*3], 3);
        fortran_read_float(fin, halodata.sp[i], 1);
        fortran_skip(fin);// time step
        fortran_read_float(fin, halodata.vir[i], 4);
        std::cout << halodata.vir[i] << std::endl;
        fortran_read_float(fin, halodata.profile[i], 2);
        fortran_skip(fin);// time step
        fortran_skip(fin);// time step
        fortran_skip(fin);// time step

    }
}

void allocate_data(Data& halodata, int ntot){
    halodata.np = new int [ntot];
    halodata.hnu = new int [ntot];
    halodata.hhost = new int [ntot * 5];
    halodata.ang = new float [ntot * 3];
    halodata.pos = new float [ntot * 3];
    halodata.vel = new float [ntot * 3];
    halodata.mass = new float [ntot];
    halodata.sp = new float [ntot];
    halodata.radius = new float [ntot * 4];
    halodata.energy = new float [ntot * 3];
    halodata.vir = new float [ntot * 4];
    halodata.profile = new float [ntot * 2];
    halodata.gal = new float [ntot * 3];
}

void load(std::string& fname, Meta& haloinfo, Data& halodata){
    
    std::fstream fhalo(fname.c_str(), std::ios::binary | std::ios::in);
    if (fhalo.fail()){
        std::cerr << "EXIT 1\n";
        std::exit(1);
    }

//    Meta haloinfo;
    load_meta(fhalo, haloinfo);
    int ntot = haloinfo.halnum + haloinfo.subnum;

//    Data2 halodata;
    allocate_data(halodata, ntot);
    load_data(fhalo, halodata, ntot);

    std::cout << halodata.pos[123] << std::endl;
    std::cout << " Load done " << std::endl;
    fhalo.close();

}


