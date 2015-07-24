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
  " usage: testDriver <GAMA_catalog>\n"
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
    if (argc < 2) {
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
    list<ComovingCoord*>::iterator i = gama_list.comovingCoordBegin();
    list<GalaxyObject*>::iterator j = gama_list.objListBegin();
    for (; i != gama_list.comovingCoordEnd(), j != gama_list.objListEnd(); ++i, ++j) {
      outf.setf(ios::fixed, ios::floatfield);
      outf << setw(11) << setprecision(6)
	   << (*j)->getRA() << " "
	   << setw(10) << setprecision(6)
	   << (*j)->getDec() << " "
	   << " "
	   << setw(5) << setprecision(3)
	   << (*j)->getRedshift() << " "
	   << "   "
	   << setw(11) << setprecision(6)
	   << (*i)->getX() << " "
	   << setw(11) << setprecision(6)
	   << (*i)->getY() << " "
	   << setw(11) << setprecision(6)
	   << (*i)->getZ() << endl;
    }


  } catch (MyException& m) {
    m.dump(cerr);
    exit(1);
  }

}
