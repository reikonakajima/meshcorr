// What does neff look like for linear power spectra?
// produce data to recreate Masahiro's fig 3.

#include "PowerSpectra.h"
#include "GrowthFunction.h"
#include "FiducialCosmology.h"
#include "AstronomicalConstants.h"
#include "Bispectra.h"
#include "SmithPS.h"
#include "LensQuantities.h"
#include "Matrix.h"
#include <iostream>
#include <iomanip>

using namespace cosmology;
using namespace std;

int
main(int argc,
     char *argv[])
{
  const double h=0.72;
  FiducialCosmology fc(0.27*h*h, 0.024, h);
  PowerLawPS lps(fc.omegaM(), fc.omegaB());
  SmithPS nlps(fc);
  HEPTBispectrum bs(nlps);

  GrowthFunction gf(fc.cosmology());
  cerr << " gr(zrec) " << gf(RecombinationRedshift) << endl;
  cerr << "Sigma8: " << lps.sigX(8/fc.H(0.))*gf(0.)/RecombinationRedshift << endl;

  double k=atof(argv[1]);
  double f = argc>2 ? atof(argv[2]) : 1.;
  double costheta = argc>3 ? atof(argv[3]) : 0.5;
  double zs = argc>4 ? atof(argv[4]) : 1.;

  if (false) {/**/
  double b, Q;
  while (cin >> k >> b >> Q) {
    cout << k
	 << " " << b
	 << " " << sqrt(bs(k,f*k,costheta,0.))*pow(k,3.)/2/PI/PI/b << endl;
  }



  cout << "DeltaSq(linear) " << lps.DeltaSq(log(k))*
    pow(gf(0.)/(1+RecombinationRedshift),2.)
       << endl;
  double dsq = nlps(k,0.)*pow(k,3.)/2/PI/PI;
  cout << "NonLinear: " << dsq << endl;
  cout << "NonLinear: " << exp(nlps.lnLinearDeltaSq(k,0.)) << endl;
  cout << "Bispectrum: " << sqrt(bs(k,f*k,costheta,0.))*pow(k,3.)/2/PI/PI << endl;
  cout << "Q: " << bs(k,f*k,costheta,0.)*
    pow(k,6.)*pow(PI,-4.)/4./dsq/dsq/3. << endl;
  exit(0);
  } /**/

  ZDistTakada zd1(1.5, 100., 0., 1.3);
  ZDistTakada zd2(1.5, 100., 1.3, 5.);
  ZDistLensWeight lw1(zd1, fc.cosmology());
  ZDistLensWeight lw2(zd2, fc.cosmology());
  //  SingleSourceLensWeight lw1(zs, fc.cosmology());
  Bkappa bk111(lw1, lw1, lw1, fc, bs);
  Bkappa bk112(lw1, lw1, lw2, fc, bs);
  Bkappa bk122(lw1, lw2, lw2, fc, bs);
  Bkappa bk222(lw2, lw2, lw2, fc, bs);

  // Set up a finite-step integration
  double dlna = 0.1;
  int nsteps = floor(log(1+zs)/dlna-0.5);
  dlna = log(1+zs) / (nsteps+0.5);
  DVector wt111(nsteps);
  DVector wt112(nsteps);
  DVector wt122(nsteps);
  DVector wt222(nsteps);
  DVector zl(nsteps+1);
  DVector davec(nsteps+1);
  for (int i=0; i<=nsteps; i++) {
    double z=exp((i+0.5)*dlna) -1.;
    davec[i] = fc.cDA(z);
    zl[i] = z;
  }
  for (int i=0; i<nsteps; i++) {
    double z=zl[i];
    double deltachi = i>0 ?
      (davec[i+1]-davec[i-1])/2. : (davec[i+1]+davec[i])/2. ;
    double deltaz = i>0 ?
      (zl[i+1]-zl[i-1])/2. : (zl[i+1]+zl[i])/2. ;
    double w1=lw1(z);
    double w2=lw2(z);
    wt111[i] = deltachi * w1*w1*w1*pow(1+z, 3.)/ fc.cDA(z);
    wt112[i] = deltachi * w1*w1*w2*pow(1+z, 3.)/ fc.cDA(z);
    wt122[i] = deltachi * w1*w2*w2*pow(1+z, 3.)/ fc.cDA(z);
    wt222[i] = deltachi * w2*w2*w2*pow(1+z, 3.)/ fc.cDA(z);
  }
  for (int i=0; i<=nsteps; i++)
    davec[i] *= HubbleLengthMpc;
  
  const double norm = pow(3*fc.omegaM()/2.,3.) * pow(HubbleLengthMpc, -6.);

  for (double logl=1.; logl<3.6; logl+=0.1) {
    double l=pow(10., logl);
    double bbin111=0.;
    double bbin112=0.;
    double bbin122=0.;
    double bbin222=0.;
    for (int i=0; i<nsteps; i++) {
      double bbs = bs(l/davec[i], f*l/davec[i], costheta, zl[i]);
      bbin111 += wt111[i] * bbs;
      bbin112 += wt112[i] * bbs;
      bbin122 += wt122[i] * bbs;
      bbin222 += wt222[i] * bbs;
    }

    bbin111 *= norm;
    bbin112 *= norm;
    bbin122 *= norm;
    bbin222 *= norm;

    const double l4=pow(l,4.)/4/PI/PI;
    cout << setw(4) << l
	 << " " << setw(8) << bk111(l, f*l, costheta)*l4
	 << " " << setw(8) << bk112(l, f*l, costheta)*l4
	 << " " << setw(8) << bk122(l, f*l, costheta)*l4
	 << " " << setw(8) << bk222(l, f*l, costheta)*l4
	 << endl;
    cout << "      "
	 << " " << setw(8) << bbin111/bk111(l, f*l, costheta)
	 << " " << setw(8) << bbin112/bk112(l, f*l, costheta)
	 << " " << setw(8) << bbin122/bk122(l, f*l, costheta)
	 << " " << setw(8) << bbin222/bk222(l, f*l, costheta)
	 << endl;
  }

}
