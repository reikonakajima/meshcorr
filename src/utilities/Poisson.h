// $Id: Poisson.h,v 2.1 2003-05-05 21:05:03 garyb Exp $
// Routines useful for Poisson distributions
#ifndef POISSON_H
#define POISSON_H

#include "Std.h"

namespace poisson {

  class PoissonError: public MyException {
  public:
    PoissonError(const string m): MyException("Poisson error; "+m)
  {}
  };

  template <class T=double>
  class Poisson {
  private:
    T mean;
  public:
    Poisson(T mean_): mean(mean_) {}
    T operator()(const int N) const;	//returns probability
    T getMean() const {return mean;}
    void setMean(const T mean_) {mean=mean_;}
    T cumulative(int N) const;		//probability of <=N
    //value of mean for which cumulative(N)=pctile
    static T percentileMean(int N, T pctile);	
  };
} // namespace poisson
#endif
