
#ifndef HISTOGRAM_H
#define HISTOGRAM_H


class HistogramError : public MyException {
 public:
 HistogramError(const string& m="") :
  MyException("HistogramError: " +m) {}
};



//
// header
//
class HistogramLogBin {
 public:
  HistogramLogBin(float rmin, float rmax, int nbin_);
  // rmin, rmax must be normalized to box side
  double operator[](int i) const { 
    if (bin) return bin[i]; 
    else throw HistogramError("HistogramLogBin not set"); }
  double rsq(int i) const {
    if (bin) return bin[i]*bin[i];
    else throw HistogramError("HistogramLogBin not set"); }
  int size() const { return Nbin; }
  double getRangeMin() const { 
    if (bin) return bin[0]; 
    else throw HistogramError("HistogramLogBin not set"); }
  double getRangeMax() const { 
    if (bin) return bin[Nbin]; 
    else throw HistogramError("HistogramLogBin not set"); }
 private:
  int    Nbin;
  double* bin;
};


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


#endif // HISTOGRAM_H
