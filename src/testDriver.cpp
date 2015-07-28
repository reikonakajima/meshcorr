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
  " usage: testDriver <GAMA_catalog> <random_catalog>\n"
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

    // read in a series of gama files
    GAMAObjectList gama_list;
    for (int i_arg=1; i_arg<argc; i_arg++) {
      const string gama_filename = argv[i_arg];
    
      /// open gama file
      ifstream gamaf(gama_filename.c_str());
      if (!gamaf) 
	throw MyException("GAMA catalog file " + gama_filename + " not found");
      gama_list.read(gamaf);
    }


    //
    // diagnostic error messages
    //
    cerr << "=== testDriver ===" << endl;
    cerr << "GAMA catalog ...... " << argv[1] << endl;
    if (argc > 2)
      for (int i=2; i<argc; ++i) {
	cerr << "             ...... " << argv[i] << endl;
      }
    cerr << "     count ........ " << gama_list.size() << endl;
    cerr << "     bounds ....... " << gama_list.getBounds() << endl;

    if (gama_list.size() == 0) {
      cerr << "no gama objects, exiting" << endl;
      return(9);
    }

    //
    // test cosmology output
    //
    double Omega_m = 0.73;
    double Omega_lambda = 0.27;
    Cosmology cosmo(Omega_m, Omega_lambda);
    gama_list.setComovingCoords(cosmo);
    string out_fname = "comovingcoord3d.out";
    cerr << "Test cosmological output in " << out_fname << endl;
    ofstream outf(out_fname.c_str());
    outf << "#ra          dec        redshift    x           y           z" << endl;
    list<GalaxyObject*>::iterator i = gama_list.objListBegin();
    for (; i != gama_list.objListEnd(); ++i) {
      outf.setf(ios::fixed, ios::floatfield);
      outf << setw(11) << setprecision(6)
	   << (*i)->getRA() << " "
	   << setw(10) << setprecision(6)
	   << (*i)->getDec() << " "
	   << " "
	   << setw(5) << setprecision(3)
	   << (*i)->getRedshift() << " "
	   << "   "
	   << setw(11) << setprecision(6)
	   << (*i)->getX() << " "
	   << setw(11) << setprecision(6)
	   << (*i)->getY() << " "
	   << setw(11) << setprecision(6)
	   << (*i)->getZ() << endl;
    }


    int nbin = 5;
    double minr = 10., maxr = 100., dx_mesh = 100.;  // Mpc/h
    HistogramLogBin rbin(minr,maxr,nbin);
    vector<double> DD, DR, RD, RR;
    vector<double> xi = LandeSzalay(gama_list, random_list, rbin, dx_mesh,
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
