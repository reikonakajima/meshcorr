// $Id: Marquardt.h,v 2.11 2002-11-05 16:44:26 garyb Exp $

// Implementation of the Marquardt-Levenberg nonlinear fitting algorithm
// in C++, and using Jarvis Matrix/Vector classes.
// The algorithm is taken from Numerical Recipes code, but has been 
// rewritten to handle more general kinds of chi-squared.

// It's a template so that classes instead of just functions can be used
// in the evaluation loops.  Not too worried about code bloat since it's
// not likely this will have multiple instantiations in one program.

// ??? Make the convergence decision a separate class and implement
// a default for it.

#ifndef MARQUARDT_H
#define MARQUARDT_H

#include "Matrix.h"
#include "Std.h"

// An exception class:
class NoConverge: public MyException {
public:
  NoConverge(const string m): MyException("Marquardt did not converge; "+m)
  {}
};

const int DefaultMaxIterations=100;
const double DefaultAbsTolerance=0.05;
const double DefaultRelTolerance=1e-4;
const double MaxLambda=1e6;

template <class P>
class MarquardtPointSum;

// Marquardt minimizer will call the object of class T with args
// (const DVector& a, double& chisq, DVector& beta, SqDMatrix& alpha))
// to build the alpha/beta/chisq from parameters a.  
// P is precision of calculation, defaults to double.
// MarquardtPointSum class is below, when you have a point function.

template <class T=MarquardtPointSum<double>, class P=double>
class Marquardt {
public: 
  explicit Marquardt(T& f):
    derivs(f), isFit(false), absTol(DefaultAbsTolerance),
    relTol(DefaultRelTolerance), bestAlpha(0),
    bestBeta(0), alpha(0), beta(0) {}
  ~Marquardt() {
    cleanup(); if (bestAlpha) delete bestAlpha;
}

  // Does the fit starting at a, returns chisq
P fit(Vector<P>& a, 
      int maxIter=DefaultMaxIterations);

  //Return (pointer to) inverse covariance matrix at last fit
  const SqMatrix<P>* getAlpha() const {return bestAlpha;}
  void setAbsTolerance(P t) {absTol=t;}
  void setRelTolerance(P t) {relTol=t;}
  P getAbsTolerance() const {return absTol;}
  P getRelTolerance() const {return relTol;}

private:
  // private copy/assignment to keep it from happening
  Marquardt(const Marquardt& m) {}
  void operator=(const Marquardt& rhs) {}
  // The object returning alpha/beta:
  T& derivs;

  // Convergence criteria:
  P absTol;
  P relTol;
  int    maxIter;

  // variables used throughout
  Vector<P> bestA;
  SqMatrix<P> *bestAlpha;		//we're going to save our best alpha
  P  bestChisq;
  P chisq;

  Vector<P> *beta;		// These are used temporarily during fit()
  Vector<P> *bestBeta;
  SqMatrix<P> *alpha;

  bool isFit;		//set flag if bestA is the best fit.
  bool lastDropWasSmall;//flag used to find 2 small drops in a row.
  void cleanup();
  bool isConverged();
};

template <class T, class P>
P
Marquardt<T,P>::fit(Vector<P>& a, 
		    int maxIter) {
  P lambda;
  int nparam = a.size();

  if (bestAlpha) delete bestAlpha;	//get proper sized matrix
  bestAlpha = new SqMatrix<P>(nparam);

  bestBeta = new DVector(nparam);
  beta = new DVector(nparam);
  alpha = new SqMatrix<P>(nparam);

  derivs(a,bestChisq,*bestBeta,*bestAlpha);//build all the matrices

  /*/ begin DEBUG
  cerr << "chisq: " << bestChisq << endl;
  cerr << "alpha: " << endl;
  cerr << *bestAlpha << endl;
  throw MyException("stop");
  // end DEBUG
  /*/

  lambda = 0.001;
  bestA = a;
  isFit = false;
  lastDropWasSmall=false;

  // If maxIter=0, the client just wanted to calculate the alpha.
  if (maxIter==0) {
    cleanup();
    return bestChisq;
  }

  // Iterate to convergence:
  for (int i=0; i<maxIter; i++) {
    // Calculate the differential step
    *alpha = *bestAlpha;
    *beta =  *bestBeta;
    for (int j=0; j<a.size(); j++)
      (*alpha)(j,j) *= (1+lambda);
    *beta /= *alpha;
    /*#ifdef DEBUG
    if (i==0) {
      cerr << "Param  Value  1st Step" << endl;
      for (int ii=0; ii<a.size(); ii++)
      cerr << setw(3) << ii 
	   << " " << setw(8) << bestA[ii] 
	   << " " << setw(8) << (*beta)[ii] 
	   << endl;
    }
    #endif*/

    a = bestA + *beta;

#ifdef DEBUG
    cerr << "a: ";
    for (int n = 0; n < a.size(); ++n)
      cerr << a[n] << " ";
    cerr << endl;
#endif

    // Get chisq and derivs at the new spot
    derivs(a,chisq,*beta,*alpha);

#ifdef DEBUG
    cerr << "Marquardt Trial " << i 
	 << " lambda=" << lambda 
	 << " chisq=" << chisq 
	 << " best chisq=" << bestChisq 
	 << endl;
#endif

    // See if we are done
    if ( isConverged() ) {
      SWAP(alpha, bestAlpha);
      cleanup();
      bestChisq = chisq;
      bestA = a;
      isFit = true;
      return bestChisq;
    }

    // Now decide circumstances of next iteration
    if (chisq < bestChisq) {
      // An improvement:
      bestChisq = chisq;
      bestA = a;
      SWAP(alpha, bestAlpha);
      SWAP(beta, bestBeta);
      lambda *= 0.1;	//Move toward quadratic solution
      // Too many times wandering down quadratic solutions, do
      // a restart with more steep descent.
      if (lambda<1e-9) lambda=1.;
    } else {
      lambda *=10.;	//Move toward steepest descent
      if (lambda > MaxLambda) {
	//cleanup();
	//a = bestA;
	//throw NoConverge("lambda is too big");
	cerr << "WARNING: Marquardt lambda=" << lambda << endl;
	SWAP(alpha, bestAlpha);
	bestChisq = chisq;
	cleanup();
	bestA = a;
	isFit = true;
	return bestChisq;
      }
    }
  }
  // Get here after max iterations
  cleanup();
  a = bestA;
  throw NoConverge("Too many iterations");
}

template <class T, class P>
void
Marquardt<T,P>::cleanup() {
      if (beta) delete beta; beta=0;
      if (alpha) delete alpha; alpha=0;
      if (bestBeta) delete bestBeta; bestBeta=0;
}

template <class T, class P>
bool
Marquardt<T,P>::isConverged() {
  P diff = bestChisq-chisq;
  if (diff<0.) return false;	//chisq must not increase
  bool isSmall = (diff < absTol) 
    || (bestChisq>0. && diff/bestChisq<relTol);
  if (isSmall && lastDropWasSmall) return true;  //two small drops=converge
  lastDropWasSmall=isSmall;
  return false;
}

// This class does what mrqmin() does in NR, namely sums of derivs, etc.,
// over a set of (x,y,sigma) data points. Precision is template arg.
// Create this object with the data and pointer to the derivative function.
template <class P=double>
class MarquardtPointSum {
 private:
  Vector<P> isig;
  const Vector<P>& x;
  const Vector<P>& y;
  const Vector<P>& sig;
  void (*func)(const Vector<P>& a,
	       const P xpt,
	       P& ypt,
	       Vector<P>& dpt);
 public:
  MarquardtPointSum(const Vector<P>& _x, 
		    const Vector<P>& _y, 
		    const Vector<P>& _sig,
		    void (*f)(const Vector<P>& a,
			      const P xpt,
			      P& ypt,
			      Vector<P>& dpt)):
    x(_x), y(_y), sig(_sig), func(f) {
    if (y.size() < x.size() || sig.size() < x.size())
      throw MyException("Data or Sigma vector size mismatch in "
			"MarquardtPointSum");
    isig.Resize(x.size());
    for (int i=0; i<x.size(); i++) isig[i] = pow(sig[i],-2.);
  }
  void operator() (const Vector<P>& a, 
		   P& chisq, 
		   Vector<P>& beta, 
		   SqMatrix<P>& alpha) const {
    chisq = 0;
    for (int i=0; i<a.size(); i++) {
      beta[i] = 0.;
      for (int j=0; j<a.size(); j++) 
	alpha(i,j)=0.;
    }
    Vector<P> derivs(a.size());
    P yfit;
    for (int k=0; k<x.size(); k++) {
      func(a, x[k], yfit, derivs);
      yfit = y[k] - yfit;
      chisq += yfit*yfit*isig[k];
      for (int i=0; i<a.size(); i++) {
	beta[i] += yfit*derivs[i]*isig[k];
	for (int j=0; j<a.size(); j++) 
	  alpha(i,j)+=derivs[i]*derivs[j]*isig[k];
      }
    }
    
    return;
  } // end operator()
};


// This is similar to MarquardtPointSum, but integrates func over the range
// (xpt-dx/2) to (xpt+dx/2) to get the ypt count.
template <class P=double>
class MarquardtHistogramSum {
 private:
  Vector<P> isig;
  const Vector<P>& x;
  const Vector<P>& y;
  const Vector<P>& sig;
  const P dx;
  void (*func)(const Vector<P>& a,
	       const P xpt,
	       P& ypt,
	       Vector<P>& dpt,
	       const P dx);
public:
  MarquardtHistogramSum(const Vector<P>& _x, 
			const Vector<P>& _y, 
			const Vector<P>& _sig,
			void (*f)(const Vector<P>& a,
				  const P xpt,
				  P& ypt,
				  Vector<P>& dpt,
				  const P dx),
			const P _dx):
    x(_x), y(_y), sig(_sig), func(f), dx(_dx) {
    if (y.size() < x.size() || sig.size() < x.size())
      throw MyException("Data or Sigma vector size mismatch in "
			"MarquardtHistogramSum");
    isig.Resize(x.size());
    for (int i=0; i<x.size(); i++) isig[i] = pow(sig[i],-2.);
  }
  void operator() (const Vector<P>& a, 
		   P& chisq, 
		   Vector<P>& beta, 
		   SqMatrix<P>& alpha) const {
    chisq = 0;
    for (int i=0; i<a.size(); i++) {
      beta[i] = 0.;
      for (int j=0; j<a.size(); j++) 
	alpha(i,j)=0.;
    }
    Vector<P> derivs(a.size());
    P yfit;
    for (int k=0; k<x.size(); k++) {
      func(a, x[k], yfit, derivs, dx);
      yfit = y[k] - yfit;
      chisq += yfit*yfit*isig[k];
      for (int i=0; i<a.size(); i++) {
	beta[i] += yfit*derivs[i]*isig[k];
	for (int j=0; j<a.size(); j++) 
	  alpha(i,j)+=derivs[i]*derivs[j]*isig[k];
      }
    }
    
    return;
  } // end operator()
};


#endif   //MARQUARDT_H
