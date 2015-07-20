#include <iostream>
#include "Std.h"
#include "NBodySimulationData.h"
using namespace nbody;
using std::setw;

const string usage =
"AutocorrelateGG: calculate galaxy autocorrelation, xi(r), given x,y,z info\n"
"                 assuming periodic boundary condition\n"
"usage:  AutocorrelateGG <x,y,z file> <ngal> <boxside> <rmin> <rmax> <nbin>\n"
"ngal:   approximate galaxy count included in file: (ngal >= actual count)\n"
"boxside, rmin, rmax:  single side, maximum/minimum r, in units of Mpc/h.\n"
"stdout: r(Mpc/h)  xi(r)\n";


int
main (int argc, char* argv[]) {


 try {

   // process arguments
   if (argc != 7) {
     cerr << usage;
     exit(1);
   }
   ifstream galf(argv[1]);
   if (!galf) {
     cerr << "file " << argv[1] << "not found" << endl;
     exit(1);
   }
   const int    objsize = atoi(argv[2]);
   const double boxside = atof(argv[3]);
   const double rmin = atof(argv[4]) / boxside;
   const double rmax = atof(argv[5]) / boxside;
   const int    nbin = atoi(argv[6]);

   // set logarithmically separated correlation bins
   HistogramLogBin rbin(rmin, rmax, nbin);

   // read in position info
   GalaxyList gallist(galf, objsize);
   /*/ DEBUG
   ofstream ofs("testout.dat");
   gallist.printPositions(ofs);
   /*/

   // setup mesh
   const int Nm = 100;   // mesh size
   GalaxyChainingMesh mesh(gallist, Nm);
   /*/ DEBUG
   ofstream ofs("testout.dat");
   mesh.printCounts(ofs);
   /*/

   // calculate autocorrelation
   double results[rbin.size()];
   mesh.calculateAutoCorrelation(rbin, gallist, results);
   // output
   for (int i = 0; i < rbin.size(); ++i) {
     cout << setw(12) << rbin[i] << " ";
     cout << setw(12) << rbin[i+1] << " ";
     cout << setw(12) << results[i] << endl;
   }
   

  } catch (MyException& m) {
    m.dump(cerr);
    exit(1);
  }

  return 0;
}

