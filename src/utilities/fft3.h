// 	$Id: fft3.h,v 1.1 2007-03-14 16:16:22 reiko3 Exp $	 
// Objects that make use of 2d FFT's in VERSION 3 of the FFTW package.
// I am ASSUMING that the FFTW is set up to do double-precision.

// Other conventions:
//	All tables have even dimensions.
//	The complex arrays (kTables) must be Hermitian, so transforms
// are real.
//	all arrays are 0-indexed.
//	rapidly varying index is 2nd one, which is y direction.
//	x-space arrays have the origin at [N/2,N/2].  This means that k
// values need to be multiplied by -1^(i+j) before/after transforms.
//	k-space arrays have origin at [0,0], since only half is stored
// anyway.  Note this means negative frequencies are at end.
//	"forward" transform, x->k, has -1 in exponent.
//	value in the table must be multiplied by "scaleby" double to get
// the correctly dimensioned/scaled value.


#ifndef FFT_H
#define FFT_H

#include "Std.h"
#include "fftw3.h"

// Class for errors
class FFTError: public MyException {
 public:
  FFTError(const string &m=""): MyException("FFT: " + m) {}
};
class FFTOutofRange: public FFTError {
 public:
  FFTOutofRange(const string &m=""): FFTError("value out of range" + m) {}
};
class FFTMalloc: public FFTError {
 public:
  FFTMalloc(const string &m=""): FFTError("malloc failure" + m) {}
};
class FFTInvalid: public FFTError {
 public:
  FFTInvalid(const string &m=""): FFTError("invalid plan or data" + m) {}
};

// FFTW3 now states that C++ complex<double> will be bit-compatible with 
// the fftw_complex type.  So all interfaces will be through our DComplex.
// And the fftw real type is now just double.

class xTable;

// Class holding k-space representation of real function.
// It will be based on an assumed Hermitian 2d square array.
// Table will be forced to be of even size.

class kTable {
  friend class xTable;
 public:
  kTable(int _N, double _dk, DComplex _value=DComplex(0.,0.));
  kTable(const kTable& rhs): array(0) {
    copy_array(rhs);
    N=rhs.N; dk=rhs.dk; valid=rhs.valid; scaleby=rhs.scaleby;};
  kTable& operator=(const kTable& rhs) {
    if (&rhs==this) return *this;
    copy_array(rhs);
    N=rhs.N; dk=rhs.dk; valid=rhs.valid; scaleby=rhs.scaleby;
    return *this;};
  ~kTable() {kill_array();};

  double getN() const {return N;}
  double getDk() const {return dk;}

  DComplex kval(int i, int j) const ;
  DComplex kval(double kx, double ky) const; //interpolate
  void kSet(int i, int j, DComplex value);
  void kSet(double kx, double ky, DComplex value); //nearest

  double xval(double x, double y) const ;

  // Fill table from a function or function object:
  void fill( DComplex func(const double kx, const double ky)) ;
  template <class T> void fill( const T &f) ;
  // New table is function of this one:
  kTable* Function( DComplex func(const double kx, 
				  const double ky, 
				  const DComplex val)) const ;
  // Integrate a function over d^2k:
  DComplex  Integrate( DComplex func(const double kx, 
				     const double ky, 
				     const DComplex val)) const ;
  // Integrate kTable over d^2k (sum of all pixels * dk * dk)
  DComplex integratePixels() const;
  // FT to give pointer to a new xTable
  xTable* Transform() const;
  void Transform(xTable& xt) const;
  // Have FFTW develop "wisdom" on doing this kind of transform
  void fftwMeasure() const;

  // Translate PSF to have origin at (x0,y0)
  void Translate(double x0, double y0);

  void operator*=(const double d) {scaleby *= d;}
  /* other things to want:  set = some constant
     translate
     convolve = multiply (with scale?)
     +=, -=, etc., with scalar or kTable.
  */
 private:
  double  scaleby;	//multiply table by this to get values
  double  dk;			//k-space increment
  int     N;			//Size in each dimension.
  bool	  valid;		//set if contains valid data right now.
  DComplex *array;		//hold the values.

  void get_array(const DComplex value=DComplex(0.,0.));	//allocate an array
  void copy_array(const kTable &rhs);	//copy an array
  void kill_array();			//deallocate array
  size_t  index(int i, int j) const;	//Return index into data array.
  // this is also responsible for bounds checking.
  void validCheck() const {if (!valid) throw FFTInvalid();};

};

// The x-space lookup table is a simple real matrix.  Force N even again,
// put origin at (N/2, N/2).
class xTable {
  friend class kTable;
 public:
  xTable(int _N, double _dx, double _value=0.);
  xTable(const xTable& rhs) {
    copy_array(rhs);
    N=rhs.N; dx=rhs.dx; valid=rhs.valid; scaleby=rhs.scaleby;};
  xTable& operator=(const xTable& rhs) {
    if (&rhs==this) return *this;
    copy_array(rhs);
    N=rhs.N; dx=rhs.dx; valid=rhs.valid; scaleby=rhs.scaleby;
    return *this;};
  ~xTable() {kill_array();};

  double xval(int i, int j) const ;
  double xval(double x, double y) const; //interpolate
  void xSet(int i, int j, double value);
  void xSet(double x, double y, double value); //nearest

  double getN() const {return N;}
  double getDx() const {return dx;}

  DComplex kval(double kx, double ky) const ;

  // Fill table from a function:
  void fill( double func(const double x, const double y)) ;
  // Integrate a (real) function over d^2x; set flag for sum:
  double  Integrate( double func(const double kx, 
				 const double ky, 
				 const double val),
		     bool sumonly=false) const ;
  // Integrate xTable over d^2x (sum of all pixels * dx * dx)
  double integratePixels() const;
  // FT to give pointer to a new kTable
  kTable* Transform() const;
  // Have FFTW develop "wisdom" on doing this kind of transform
  void fftwMeasure() const;

  /* other things to want:  set = some constant
     translate
     transform to kTable
     convolve (via kTable?)
     +=, -=, etc., with scalar or kTable.
  */
 private:
  double  scaleby;	//multiply table by this to get values
  double  dx;			//k-space increment
  int     N;			//Size in each dimension.
  bool	  valid;		//set if contains valid data right now.
  double  *array;		//hold the values.

  void get_array(const double value);	//allocate an array
  void copy_array(const xTable &rhs);	//copy an array
  void kill_array();			//deallocate array
  size_t  index(int i, int j) const;	//Return index into data array.
  // this is also responsible for bounds checking.
  void validCheck() const {if (!valid) throw FFTInvalid();};
};

// Fill table from a function class:
template <class T>
void
kTable::fill( const T& f) {
  DComplex *zptr=array;
  double kx, ky;
  for (int i=0; i< N/2; i++) {
   ky = i*dk;
    for (int j=0; j< N/2+1 ; j++) {
      kx = j*dk;
      *(zptr++) = f(kx,ky);
    }
  }
  // wrap to the negative ky's
  for (int i=-N/2; i< 0; i++) {
   ky = i*dk;
    for (int j=0; j< N/2+1 ; j++) {
      kx = j*dk;
      *(zptr++) = f(kx,ky);
    }
  }
  return;
}

#endif
