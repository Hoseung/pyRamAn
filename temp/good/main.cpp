#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include "c_rd_halo.h"

int main(int argc, char *argv[])
{
    Meta haloinfo;
    Data2 halodata;

//    std::string fname ("/home/hoseung/Work/data/29176/GalaxyMaker/gal/tree_bricks187");
    std::string fname ("/home/hoseung/Work/data/29176/halo/DM/tree_bricks187");
    load(fname, haloinfo, halodata);
    std::cout << " Done " << std::endl;
}
