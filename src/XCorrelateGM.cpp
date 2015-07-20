#include <iostream>
#include "Std.h"
#include "NBodySimulationData.h"
using namespace nbody;
using std::setw;
using std::ios_base;

const string usage =
"XCorrelateGM: calculate galaxy-matter cross correlation, xi_gm(r), \n"
"              given x,y,z info, assuming periodic boundary condition\n"
"usage:  XCorrelate <dmfile> <galfile> <ngal> <boxside> <rmin> <rmax> <nbin>\n"
"ngal:   approximate galaxy count included in file: (ngal >= actual count)\n"
"boxside, rmin, rmax:  single side, maximum/minimum r, in units of [Mpc/h]\n"
"stdout: r(Mpc/h)  xi(r)\n";


int
main (int argc, char* argv[]) {


 try {

   // process arguments
   if (argc != 8) {
     cerr << usage;
     exit(1);
   }
   int i = 0;
   ifstream dmf(argv[++i]);
   if (!dmf) {
     cerr << "file " << argv[i] << " not found" << endl;
     exit(1);
   }
   dmf.close();
   ifstream galf(argv[++i]);
   if (!galf) {
     cerr << "file " << argv[i] << " not found" << endl;
     exit(1);
   }
   const int    objsize = atoi(argv[++i]);
   const double boxside = atof(argv[++i]);
   const double rmin = atof(argv[++i]) / boxside;
   const double rmax = atof(argv[++i]) / boxside;
   const int    nbin = atoi(argv[++i]);

   // set logarithmically separated correlation bins
   HistogramLogBin rbin(rmin, rmax, nbin);

   // read in position info
   DMParticleList dmlist(argv[1]);
   //GalaxyList gallist(galf, objsize);

   // setup mesh
   const int Nm = 100;   // mesh size
   //GalaxyChainingMesh mesh(gallist, Nm);
   /*/ DEBUG
   ofstream ofs("testout.dat");
   mesh.printCounts(ofs);
   /*/

   /*/ calculate autocorrelation
   double results[rbin.size()];
   mesh.calculateAutoCorrelation(rbin, gallist, results);
   // output
   for (int i = 0; i < rbin.size(); ++i) {
     cout << setw(12) << rbin[i] << " ";
     cout << setw(12) << rbin[i+1] << " ";
     cout << setw(12) << results[i] << endl;
   }
   /*/

  } catch (MyException& m) {
    m.dump(cerr);
    exit(1);
  }

  std::cerr << "Strawberry!\n";
  return 0;
}

