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
  " usage: testDriver <smass_bin_name> <min_log10_smass> <max_log10_smass>\n"
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
    if (argc < 4) {
      cerr << usage;
      exit(2);
    }

    int i_arg = 0;
    const string smass_bin_name = argv[++i_arg];
    const float min_log10_smass = atof(argv[++i_arg]);
    const float max_log10_smass = atof(argv[++i_arg]);

    /// read in the 3 gama files (file specific to the Bonn Euclid machines)
    vector<string> gama_fits_filename(3);
    gama_fits_filename[0] = "~/1project/galbias/wp/data_lens_3d/G09_3d.fits";
    gama_fits_filename[1] = "~/1project/galbias/wp/data_lens_3d/G12_3d.fits";
    gama_fits_filename[2] = "~/1project/galbias/wp/data_lens_3d/G15_3d.fits";
    GAMAObjectList master_gama_list;
    for (int i=0; i<gama_fits_filename.size(); ++i) {
      master_gama_list.read(gama_fits_filename[i]);
    }
    cerr << "GAMA data read" << endl;
    GAMAObjectList gama_list = master_gama_list.cullByLogMStar(min_log10_smass, max_log10_smass);

    /// read in the 3 random file
    vector<string> random_fits_filename(3);
    random_fits_filename[0] = "~/1project/galbias/wp/rand_30x/G09.rand.fits";
    random_fits_filename[1] = "~/1project/galbias/wp/rand_30x/G12.rand.fits";
    random_fits_filename[2] = "~/1project/galbias/wp/rand_30x/G15.rand.fits";

    GAMAObjectList master_random_list;
    for (int i=0; i<random_fits_filename.size(); ++i) {
      master_random_list.read(random_fits_filename[i]);
    }
    cerr << "master random data read ";
    int decimate_factor = 10;
    GalaxyObjectList random_list = master_random_list.cull(decimate_factor);
    cerr << "and culled" << endl;


    //
    // diagnostic error messages
    //
    cerr << "=== testDriver ===" << endl;
    for (int i=0; i<gama_fits_filename.size(); ++i) {
      cerr << "GAMA catalog ...... " << gama_fits_filename[i] << endl;
    }
    cerr << "     count ........ " << gama_list.size() << endl;
    cerr << "     bounds ....... " << gama_list.getBounds() << endl;
    for (int i=0; i<random_fits_filename.size(); ++i) {
      cerr << "Random catalog .... " << random_fits_filename[i] << endl;
    }
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
    string out_filename = smass_bin_name+".xi.out";
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
