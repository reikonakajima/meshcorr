
#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include "Std.h"


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


#endif // HISTOGRAM_H
