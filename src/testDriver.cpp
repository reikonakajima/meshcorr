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
    int n = 10;
    bool periodic = false;
    vector<GalaxyObject*> gamavector = gama_list.getVectorForm();
    double xmin, xmax, ymin, ymax, zmin, zmax;
    bool addEpsilon = false;
    gama_list.getXYZMinMax(xmin, xmax, ymin, ymax, zmin, zmax, addEpsilon);
    cerr << xmin << " " << xmax << " "
	 << ymin << " " << ymax << " "
	 << zmin << " " << zmax << endl;
    cerr << static_cast<int>(2.5);
    Mesh<GalaxyObject*, double> mesh(n, n, n, gamavector, periodic,
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
    double xsize = (xmax - xmin);
    double ysize = (ymax - ymin);
    double zsize = (zmax - zmin);
    cerr << xsize << " " << ysize << " " << zsize << endl;
    int ix0 = static_cast<int>(n / xsize * (x0-xmin));
    int iy0 = static_cast<int>(n / ysize * (y0-ymin));
    int iz0 = static_cast<int>(n / zsize * (z0-zmin));
    cerr << ix0 << " " << iy0 << " " << iz0 << endl;
    // find maximum search radius in the x,y,z axis (assume cartesian coordinates)
    int srx = static_cast<int>(maxr / xsize) + 2;
    int sry = static_cast<int>(maxr / ysize) + 2;
    int srz = static_cast<int>(maxr / zsize) + 2;
    cerr << srx << " " << sry << " " << srz << endl;
    int srx_min = (ix0-srx) < 0 ? 0 : (ix0-srx);
    int srx_max = (ix0+srx) >= n ? n-1 : (ix0+srx);
    int sry_min = (iy0-sry) < 0 ? 0 : (iy0-sry);
    int sry_max = (iy0+sry) >= n ? n-1 : (iy0+sry);
    int srz_min = (iz0-srz < 0) ? 0 : (iz0-srz);
    int srz_max = (iz0+srz >= n) ? n-1 : (iz0+srz);
    cerr << srx_min << " " << sry_min << " " << srz_min << endl;
    cerr << srx_max << " " << sry_max << " " << srz_max << endl;
    // identify all mesh within search radius
    for (int ix = srx_min; ix <= srx_max; ix++) {
      for (int iy = sry_min; iy <= sry_max; iy++) {
	for (int iz = srz_min; iz <= srz_max; iz++) {
	  int ii;
	  if (ix*ix + iy*iy + iz*iz < srx) {
	    ;
	  }
	}
      }
    }

      // calculate mesh id

      // calculate distance to all galaxies

      // include in neighbor list if within range

  } catch (MyException& m) {
    m.dump(cerr);
    exit(1);
  }

}
