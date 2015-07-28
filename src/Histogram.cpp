#include "Histogram.h"

//
// implementation
//
HistogramLogBin::HistogramLogBin(float rmin, float rmax, int nbin_) 
: Nbin(nbin_) {
  if (rmin <= 0. || rmax <= 0.) 
    throw HistogramError("negative or zero rmin or rmax");
  /*  // only needed for periodic conditions
  if (rmin > 1. || rmax > 1.)
    throw HistogramError("rmin or rmax not normalized to box size");
  */

  bin = new double[Nbin+1];
  double logoffset = log(rmin);
  double logsize = log(rmax/rmin);
  double dlogsize = logsize / Nbin;
  
  bin[0] = rmin;
  for (int i = 1; i < (Nbin+1); ++i) {
    bin[i] = exp( logoffset + i*dlogsize );
  }
}


