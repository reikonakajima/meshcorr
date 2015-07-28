//
// meshDriver.cpp   : For testing the Mesh C++ module
//
#include <iostream>
#include <iomanip>
#include "Std.h"
#include "StringStuff.h"
#include "Bounds.h"
#include "GAMAObjects.h"
#include "Cosmology.h"
#include "Mesh.h"
#include "Histogram.h"
#include <CCfits/CCfits>
using std::ios;
using std::ostringstream;
using std::setw;
using std::setfill;
using std::fixed;
using std::setprecision;


const string usage =
  "\n"
  "testDriver: calculate tangential shear around a point (galaxy-galaxy lensing signal)\n"
  "\n"
  " usage: testDriver <GAMA_FITS_catalog> <random_FITS_catalog>\n"
  "  GAMA_catalog:  lens catalog\n"
  "  \n"
  //  " output #1: file name:\" "+outfprefix+suffix+"\"\n"
  " stdin:  (none)\n"
  " stdout: (none)\n";


int
main(int argc, char* argv[]) {

  try {

    //
    // process arguments
    //
    if (argc < 3) {
      cerr << usage;
      exit(2);
    }

    // read in gama files
    int i_arg = 0;
    const string gama_filename = argv[++i_arg];
    ifstream gamaf(gama_filename.c_str());
    if (!gamaf)
      throw MyException("GAMA catalog file " + gama_filename + " not found");
    GAMAObjectList gama_list(gama_filename);

    /*
    /// open random file
    const string random_filename = argv[++i_arg];
    ifstream randomf(gama_filename.c_str());
    if (!randomf)
      throw MyException("Random catalog file " + random_filename + " not found");
    GalaxyObjectList random_list(randomf);
    */

    //
    // diagnostic error messages
    //
    cerr << "=== testDriver ===" << endl;
    cerr << "GAMA catalog ...... " << argv[1] << endl;
    cerr << "     count ........ " << gama_list.size() << endl;
    cerr << "     bounds ....... " << gama_list.getBounds() << endl;
    cerr << "Random catalog .... " << argv[2] << endl;
    /*
    cerr << "     count ........ " << random_list.size() << endl;
    cerr << "     bounds ....... " << random_list.getBounds() << endl;
    */

    if (gama_list.size() == 0) {
      cerr << "no gama objects, exiting" << endl;
      return(9);
    }

    //
    // set XYZ accoring to cosmology
    //
    double Omega_m = 0.73;
    double Omega_lambda = 0.27;
    Cosmology cosmo(Omega_m, Omega_lambda);
    cerr << "calculating comoving coordinates for data" << endl;
    gama_list.setComovingCoords(cosmo);
    //cerr << "calculating comoving coordinates for randoms" << endl;
    // random_list.setComovingCoords(cosmo);



    int nbin = 5;
    double minr = 10., maxr = 100., dx_mesh = 100.;  // Mpc/h
    HistogramLogBin rbin(minr,maxr,nbin);
    vector<double> DD, DR, RD, RR;
    vector<double> xi = LandeSzalay(gama_list, gama_list, rbin, dx_mesh,
				    DD, DR, RD, RR);

    for (int i = 0; i < rbin.size(); ++i) {
      // DEBUG OUTPUT
      cerr << i << " : "
	   << rbin[i] << "," << rbin[i+1] << " : "
	   << DD[i] << " "
	   << DR[i] << " "
	   << RD[i] << " "
	   << RR[i] << " "
	   << xi[i] << endl;
    }


  } catch (MyException& m) {
    m.dump(cerr);
    exit(1);
  }

}
