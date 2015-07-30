//
// gama_smass_bin_3dcorr.cpp   : For testing the Mesh C++ module
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
using std::scientific;
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
    if (argc < 2) {
      cerr << usage;
      exit(2);
    }

    const string data_dir = argv[1];

    /// read in the 3 gama files (only runs on the Bonn Euclid machines)
    string gama_fits_filename[3];
    gama_fits_filename[0] = "~/1project/galbias/wp/data_lens_3d/G09_3d.fits";
    gama_fits_filename[1] = "~/1project/galbias/wp/data_lens_3d/G12_3d.fits";
    gama_fits_filename[2] = "~/1project/galbias/wp/data_lens_3d/G15_3d.fits";

    GAMAObjectList gama_list;
    cerr << "GAMA data read" << endl;


    /// read in the 3 random file
    const string random_fits_filename = argv[++i_arg];
    ifstream randomf(random_fits_filename.c_str());
    if (!randomf)
      throw MyException("Random catalog file " + random_fits_filename + " not found");
    GAMAObjectList master_random_list(random_fits_filename);
    cerr << "master random data read ";
    int decimate_factor = 5;
    GalaxyObjectList random_list = master_random_list.cull(decimate_factor);
    cerr << "and culled" << endl;


    //
    // diagnostic error messages
    //
    cerr << "=== testDriver ===" << endl;
    cerr << "GAMA catalog ...... " << argv[1] << endl;
    cerr << "     count ........ " << gama_list.size() << endl;
    cerr << "     bounds ....... " << gama_list.getBounds() << endl;
    cerr << "Random catalog .... " << argv[2] << endl;
    cerr << "     count ........ " << random_list.size() << endl;
    cerr << "     bounds ....... " << random_list.getBounds() << endl;

    if (gama_list.size() == 0) {
      cerr << "no gama objects, exiting" << endl;
      return(9);
    }

    if (random_list.size() == 0) {
      cerr << "no random objects, exiting" << endl;
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
    cerr << "calculating comoving coordinates for randoms" << endl;
    random_list.setComovingCoords(cosmo);



    int nbin = 40;
    double minr = 0.1 , maxr = 100., dx_mesh = 20.;  // Mpc/h
    HistogramLogBin rbin(minr,maxr,nbin);
    vector<double> DD, DR, RD, RR, mean_r;
    bool isAutoCorr = true;
    vector<double> xi = LandeSzalay(gama_list, random_list, rbin, dx_mesh,
				    DD, DR, RD, RR, mean_r, isAutoCorr);

    // save to output file
    string out_filename = gama_fits_filename+".xi.out";
    ofstream outf(out_filename.c_str());
    outf << "#i    rmin[i]  rmax[i]  mean_r[i]       DD             DR             RD             RR             xi[i]" << endl;
    for (int i = 0; i < rbin.size(); ++i) {
      // OUTPUT
      outf << fixed << setw(3)
	   << i << "  "
	   << fixed << setw(8) << setprecision(4)
	   << rbin[i] << " "
	   << fixed << setw(8) << setprecision(4)
	   << rbin[i+1] << "   "
	   << fixed << setw(8) << setprecision(4)
	   << mean_r[i] << "  "
	   << setw(14) << setprecision(6) << scientific
	   << DD[i] << " "
	   << setw(14) << setprecision(6)
	   << DR[i] << " "
	   << setw(14) << setprecision(6)
	   << RD[i] << " "
	   << setw(14) << setprecision(6)
	   << RR[i] << " "
	   << setw(14) << setprecision(6)
	   << xi[i] << endl;
    }


  } catch (MyException& m) {
    m.dump(cerr);
    exit(1);
  }

}
