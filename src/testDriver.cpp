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

    //
    // test Mesh compatibility
    //
    double dx = 100.;  // size of the mesh cell
    bool periodic = false;
    vector<GalaxyObject*> gamavector = gama_list.getVectorForm();
    cerr << "gama_list size: "  << gama_list.size() << endl;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    bool addEpsilon = true;
    gama_list.getXYZMinMax(xmin, xmax, ymin, ymax, zmin, zmax, addEpsilon);
    cerr << xmin << " " << xmax << " "
	 << ymin << " " << ymax << " "
	 << zmin << " " << zmax << endl;
    cerr << "dx=dy=dz: " << dx << endl;
    Mesh<GalaxyObject*, double> mesh(dx, dx, dx, gamavector, periodic,
				     xmin, xmax, ymin, ymax, zmin, zmax);


    //
    // test Mesh nearest neighbor of x0,y0,z0  (getNeighborList)
    //
    double x0 = -300.;
    double y0 = 400.;
    double z0 = 0.;
    cerr << x0 << " " << y0 << " " << z0 << endl;
    double maxr = 100.;
    double minr = 0.;
    // find mesh position for current object x0,y0,z0
    list<int> nbr_index = mesh.getNearMeshList(x0, y0, z0, maxr, minr);
    cerr << "nbr_index size: "  << nbr_index.size() << endl;

    string out_fname_test = "test.out";
    ofstream outf_test(out_fname_test.c_str());
    outf_test << "#ra          dec        redshift    x           y           z" << endl;
    for (std::list<int>::iterator ii=nbr_index.begin(); ii!=nbr_index.end(); ii++) {
      outf_test.setf(ios::fixed, ios::floatfield);
      outf_test << setw(11) << setprecision(6)
	   << gamavector[*ii]->getRA() << " "
	   << setw(10) << setprecision(6)
	   << gamavector[*ii]->getDec() << " "
	   << " "
	   << setw(5) << setprecision(3)
	   << gamavector[*ii]->getRedshift() << " "
	   << "   "
	   << setw(11) << setprecision(6)
	   << gamavector[*ii]->getX() << " "
	   << setw(11) << setprecision(6)
	   << gamavector[*ii]->getY() << " "
	   << setw(11) << setprecision(6)
	   << gamavector[*ii]->getZ() << endl;
    }

    /*
      // calculate mesh id

      // calculate distance to all galaxies

      // include in neighbor list if within range
     */

  } catch (MyException& m) {
    m.dump(cerr);
    exit(1);
  }

}
