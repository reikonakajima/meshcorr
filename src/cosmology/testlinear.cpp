// What does neff look like for linear power spectra?

#include "PowerSpectra.h"
#include "GrowthFunction.h"
#include "FiducialCosmology.h"
#include "AstronomicalConstants.h"
#include "Bispectra.h"
#include "SmithPS.h"
#include "LensQuantities.h"
#include <iostream>
#include <iomanip>

using namespace cosmology;
using namespace std;

int
main(int argc,
     char *argv[])
{
  FiducialCosmology fc;
  PowerLawPS lps(fc.omegaM(), fc.omegaB());
  SmithPS nlps(fc);
  HEPTBispectrum bs(nlps);


  double l1=atof(argv[1]);
  double f = argc>2 ? atof(argv[2]) : 1.;
  double costheta = argc>3 ? atof(argv[3]) : 0.5;
  double zs = argc>4 ? atof(argv[4]) : 1.;

  SingleSourceLensWeight lw1(zs, fc.cosmology());
  Bkappa bk(lw1, lw1, lw1, fc, bs);

  for (double logl=1.; logl<3.5; logl+=0.25) {
    double l=pow(10., logl);
    cout << setw(4) << l
	 << " " << setw(8) << bk(l, f*l, costheta)
	 << " " << setw(8) << bk(l, f*l, costheta)*pow(l,4.)/4/PI/PI
	 << endl;
  }

}
