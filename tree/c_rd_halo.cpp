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
    fin.read((char *) &dummy, sizeof(dummy));
    return ;
}

void fortran_read_long(std::fstream& fin, long &data, int len){
    int dummy;
    fin.read((char *) &dummy, sizeof(dummy));
    fin.read((char *) &data, sizeof(long) * len);
    fin.read((char *) &dummy, sizeof(dummy));
    return ;
}

void fortran_read_float(std::fstream& fin, float &data, int len){
    int dummy;
//    double data;
    fin.read((char *) &dummy, sizeof(dummy));
    fin.read((char *) &data, sizeof(float) * len);
    fin.read((char *) &dummy, sizeof(dummy));
    return ;
}

void fortran_read_double(std::fstream& fin, double &data, int len){
    int dummy;
//    double data;
    fin.read((char *) &dummy, sizeof(dummy));
    fin.read((char *) &data, sizeof(double) * len);
    fin.read((char *) &dummy, sizeof(dummy));
    return ;
}

void load_meta(std::fstream& fin, Meta& haloinfo){
    int dummy [2];
    fortran_read_int(fin, haloinfo.nbodies, 1);
    fortran_read_float(fin, haloinfo.massp, 1);
    fortran_read_float(fin, haloinfo.aexp, 1);
    fortran_read_float(fin, haloinfo.omegat, 1);
    fortran_read_float(fin, haloinfo.age, 1);
    fortran_read_int(fin, dummy[0], 2);
    haloinfo.halnum = dummy[0];
    haloinfo.subnum = dummy[1];
}

void load_meta_d(std::fstream& fin, Meta2& haloinfo){
    int dummy [2];
    fortran_read_int(fin, haloinfo.nbodies, 1);
    fortran_read_double(fin, haloinfo.massp, 1);
    fortran_read_double(fin, haloinfo.aexp, 1);
    fortran_read_double(fin, haloinfo.omegat, 1);
    fortran_read_double(fin, haloinfo.age, 1);
    fortran_read_int(fin, dummy[0], 2);
    haloinfo.halnum = dummy[0];
    haloinfo.subnum = dummy[1];
}


void load_data_gal(std::fstream& fin, Meta& haloinfo,  Data& halodata, int tothal, int is_gal, int read_mbp){
    //int* ids;
    int iskip=0;
    for (int i=0; i < tothal; i++){
        fortran_read_int(fin, halodata.np[i], 1);
        //ids = new int [halodata.np[i]];
        fortran_read_int(fin, haloinfo.allID[iskip], halodata.np[i]);
        //fortran_read_int(fin, ids[0], halodata.np[i]);
        fortran_read_int(fin, halodata.hnu[i], 1);
        fortran_skip(fin);// time step
        fortran_read_int(fin, halodata.hhost[i*5], 5);
        if(read_mbp>0)
      	{
      	  fortran_read_long(fin, halodata.imbp[i], 1);
          if (i < 100){
            std::cout<<halodata.imbp[i]<<std::endl;
          }
      	}
        fortran_read_float(fin, halodata.mass[i], 1);
        fortran_read_float(fin, halodata.pos[i*3], 3);
        fortran_read_float(fin, halodata.vel[i*3], 3);
        fortran_read_float(fin, halodata.ang[i*3], 3);
        fortran_read_float(fin, halodata.radius[i*4], 4);
        fortran_read_float(fin, halodata.energy[i*3], 3);
        fortran_read_float(fin, halodata.sp[i], 1);
        if (is_gal > 0)
        {
            fortran_read_float(fin, halodata.gal[i*3], 3);
        }
        fortran_read_float(fin, halodata.vir[i*4], 4);
        fortran_read_float(fin, halodata.profile[i*2], 2);
        if (is_gal > 0)
        {
            fortran_read_int(fin, halodata.g_nbin[i], 1);
            fortran_read_float(fin, halodata.g_rr[i*100], 100);
            fortran_read_float(fin, halodata.g_rho[i*100], 100);
        }
        iskip += halodata.np[i];
    }
}

void load_data_gal_d(std::fstream& fin, Meta2& haloinfo,  Data2& halodata, int tothal, int is_gal, int read_mbp){
    //int* ids;
    int iskip=0;
    for (int i=0; i < tothal; i++){
        fortran_read_int(fin, halodata.np[i], 1);
        //ids = new int [halodata.np[i]];
        fortran_read_int(fin, haloinfo.allID[iskip], halodata.np[i]);
        //fortran_read_int(fin, ids[0], halodata.np[i]);
        fortran_read_int(fin, halodata.hnu[i], 1);
        fortran_skip(fin);// time step
        fortran_read_int(fin, halodata.hhost[i*5], 5);
        fortran_read_double(fin, halodata.mass[i], 1);
        fortran_read_double(fin, halodata.pos[i*3], 3);
        fortran_read_double(fin, halodata.vel[i*3], 3);
        fortran_read_double(fin, halodata.ang[i*3], 3);
        fortran_read_double(fin, halodata.radius[i*4], 4);
        fortran_read_double(fin, halodata.energy[i*3], 3);
        fortran_read_double(fin, halodata.sp[i], 1);
        if (is_gal > 0)
        {
            fortran_read_double(fin, halodata.gal[i*3], 3);
        }
        fortran_read_double(fin, halodata.vir[i*4], 4);
        fortran_read_double(fin, halodata.profile[i*2], 2);
        if (is_gal > 0)
        {
            fortran_read_int(fin, halodata.g_nbin[i], 1);
            fortran_read_double(fin, halodata.g_rr[i*100], 100);
            fortran_read_double(fin, halodata.g_rho[i*100], 100);
        }
        iskip += halodata.np[i];
    }
}


void allocate_data(Data& halodata, int ntot, int is_gal, int read_mbp){
    halodata.np = new int [ntot];
    halodata.hnu = new int [ntot];
    halodata.hhost = new int [ntot * 5];
    halodata.imbp = new long [ntot];
    halodata.ang = new float [ntot * 3];
    halodata.pos = new float [ntot * 3];
    halodata.vel = new float [ntot * 3];
    halodata.mass = new float [ntot];
    halodata.sp = new float [ntot];
    halodata.radius = new float [ntot * 4];
    halodata.energy = new float [ntot * 3];
    halodata.vir = new float [ntot * 4];
    halodata.profile = new float [ntot * 2];
    if (is_gal > 0){
        halodata.gal = new float [ntot * 3];
        halodata.g_nbin = new int [ntot];
        halodata.g_rr = new float [ntot*100];  // rr and rho are 100-element arrays
        halodata.g_rho = new float [ntot*100];
    }
}

void allocate_data_d(Data2& halodata, int ntot, int is_gal, int read_mbp){
    halodata.np = new int [ntot];
    halodata.hnu = new int [ntot];
    halodata.hhost = new int [ntot * 5];
    halodata.imbp = new long [ntot];
    halodata.ang = new double [ntot * 3];
    halodata.pos = new double [ntot * 3];
    halodata.vel = new double [ntot * 3];
    halodata.mass = new double [ntot];
    halodata.sp = new double [ntot];
    halodata.radius = new double [ntot * 4];
    halodata.energy = new double [ntot * 3];
    halodata.vir = new double [ntot * 4];
    halodata.profile = new double [ntot * 2];
    if (is_gal > 0){
        halodata.gal = new double [ntot * 3];
        halodata.g_nbin = new int [ntot];
        halodata.g_rr = new double [ntot*100];  // rr and rho are 100-element arrays
        halodata.g_rho = new double [ntot*100];
    }
}


void load(std::string& fname, Meta& haloinfo, Data& halodata, int nbodies, int is_gal, int read_mbp){
    std::fstream fhalo(fname.c_str(), std::ios::binary | std::ios::in);
    if (fhalo.fail()){
        std::cerr << "EXIT 1\n";
        std::exit(1);
    }
//    Meta haloinfo;
    haloinfo.allID = new int [nbodies];
    load_meta(fhalo, haloinfo);
    int ntot = haloinfo.halnum + haloinfo.subnum;

//    Data2 halodata;
    allocate_data(halodata, ntot, is_gal, read_mbp);
    load_data_gal(fhalo, haloinfo, halodata, ntot, is_gal, read_mbp);

//    std::cout << halodata.pos[123] << std::endl;
//    std::cout << " Load done " << std::endl;
    fhalo.close();
}

void load_d(std::string& fname, Meta2& haloinfo, Data2& halodata, int nbodies, int is_gal, int read_mbp){
    std::fstream fhalo(fname.c_str(), std::ios::binary | std::ios::in);
    if (fhalo.fail()){
        std::cerr << "EXIT 1\n";
        std::exit(1);
    }
//    Meta haloinfo;
    haloinfo.allID = new int [nbodies];
    load_meta_d(fhalo, haloinfo);
    int ntot = haloinfo.halnum + haloinfo.subnum;

//    Data2 halodata;
    allocate_data_d(halodata, ntot, is_gal, read_mbp);
    load_data_gal_d(fhalo, haloinfo, halodata, ntot, is_gal, read_mbp);

//    std::cout << halodata.pos[123] << std::endl;
//    std::cout << " Load done " << std::endl;
    fhalo.close();
}
