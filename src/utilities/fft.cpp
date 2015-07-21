// Routines for FFTW interface objects
// This time to use Version 3 of FFTW.
// 	$Id: fft.cpp,v 2.6 2007-02-11 15:24:25 garyb Exp $	

#include "fft.h"

#ifndef PI
#define PI 3.1415926535
#endif

kTable::kTable(int _N, 
	       double _dk, 
	       DComplex value) {
  N = 2*((_N+1)/2);	//Round size up to even.
  dk = _dk;
  get_array(value);
  valid = true;
  scaleby=1.;
  return;
}

size_t
kTable::index(int i, int j) const {
  // adjust for xTable with origin in center.
  if (i<-N/2 || i>N/2 || j<-N/2 || j>N/2) throw FFTOutofRange() ;
  if (j<0) {
    j=-j; i=-i;	//need the conjugate in this case
  }
  if (i<0) i+=N;
  return i*(N/2+1)+j;
}

DComplex kTable::kval(int i, int j) const {
  validCheck();
  DComplex retval=scaleby*array[index(i,j)];
  if (j<0) return conj(retval);
  else return retval;
}

void
kTable::kSet(int i, int j, DComplex value) {
  validCheck();
  if (j<0) {
    array[index(i,j)]=conj(value)/scaleby;
    if (j==-N/2) array[index(-i,j)]=value/scaleby;
  }
  else {
    array[index(i,j)]=value/scaleby;
    if (j==0 || j==N/2) array[index(-i,j)]=conj(value)/scaleby;
  }
  return;
}

void
kTable::get_array(const DComplex value) {
  array = (DComplex*) fftw_malloc(sizeof(DComplex)*N*(N/2+1));
  for (int i=0; i<N*(N/2+1); i++)
    array[i]=value;
  return;
}

void
kTable::copy_array(const kTable& rhs) {
  if (array!=0 && N!=rhs.N) {
    fftw_free(array);
    array = 0;
  }
  N = rhs.N;
  if (!rhs.valid || rhs.array==0) {
    valid=false;
    fftw_free(array);
    array = 0;
    return;
  } else
    valid = true;
  if (array==0) array = (DComplex*) fftw_malloc(sizeof(DComplex)*N*(N/2+1));

  for (int i=0; i<N*(N/2+1); i++)
    array[i]=rhs.array[i];
  return;
}

void
kTable::kill_array() {
  fftw_free(array);
  array=0;
  return;
}

// Interpolate table (linearly) to some specific k:
DComplex 
kTable::kval(double kx, double ky) const {
  kx /= dk;
  ky /= dk;
  int i = static_cast<int> (floor(ky));
  double fy = ky-i;
  int j = static_cast<int> (floor(kx));
  double fx = kx-j;

  return (1-fx)*( (1-fy)*kval(i,j) + fy*kval(i+1,j)) +
    fx*( (1-fy)*kval(i,j+1) + fy*kval(i+1,j+1));
}

// Set table point nearest to some k
void
kTable::kSet(double kx, double ky, DComplex value) {
  kx /= dk;
  ky /= dk;
  int j = static_cast<int> (floor(kx+0.5));
  int i = static_cast<int> (floor(ky+0.5));
  kSet(i,j,value);
  return;
}
  
// Fill table from a function:
void
kTable::fill( DComplex func(const double kx, const double ky)) {
  DComplex *zptr=array;
  double kx, ky;
  DComplex *tmp1 = new DComplex[N/2];
  DComplex *tmp2 = new DComplex[N/2];
  /*
  for (int i=0; i< N/2; i++) {
   ky = i*dk;
    for (int j=0; j< N/2+1 ; j++) {
      kx = j*dk;
      *(zptr++) = func(kx,ky);
    }
  }
  // wrap to the negative ky's
  for (int i=-N/2; i< 0; i++) {
   ky = i*dk;
    for (int j=0; j< N/2+1 ; j++) {
      kx = j*dk;
      *(zptr++) = func(kx,ky);
    }
  }
  */

  // [ky/dk] = i = 0
  for (int j=0; j< N/2+1 ; j++) {
    kx = j*dk;
    *(zptr++) = func(kx,0);                  // [kx/dk] = j = 0 to N/2
  }
  // [ky/dk] = i = 1 to (N/2-1)
  for (int i=1; i< N/2; i++) {
    ky = i*dk;
    *(zptr++) = tmp1[i] = func(0,ky);        // [kx/dk] = j = 0
    for (int j=1; j< N/2 ; j++) {    
      kx = j*dk;
      *(zptr++) = func(kx,ky);               // [kx/dk] = j = 1 to (N/2-1)
    }
    *(zptr++) = tmp2[i] = func((N/2)*dk,ky); // [kx/dk] = j =N/2
  }
  // Wrap to the negative ky's
  // [ky/dk] = i = -N/2
  for (int j=0; j< N/2+1 ; j++) {
    kx = j*dk;
    *(zptr++) = func(kx,-N/2);         // [kx/dk] = j = 0 to N/2   
  }
  // [ky/dk] = i = (-N/2+1) to (-1)
  for (int i=-N/2+1; i< 0; i++) {
   ky = i*dk;
   *(zptr++) = conj(tmp1[-i]);       // [kx/dk] = j = 0
    for (int j=1; j< N/2 ; j++) {
      kx = j*dk;
      *(zptr++) = func(kx,ky);         // [kx/dk] = j = 1 to (N/2-1)
    }
    *(zptr++) = conj(tmp2[-i]);      // [kx/dk] = j = N/2
  }
  return;
}


// Integrate a function over k - can be function of k or of PSF(k)
DComplex
kTable::Integrate( DComplex func(const double kx, 
				 const double ky, 
				 const DComplex val)) const {
  DComplex sum=0.;
  DComplex val;
  double kx, ky;
  DComplex *zptr=array;
  // Do the positive y frequencies
  for (int i=0; i<= N/2; i++) {
    ky = i*dk;
    val = *(zptr++) * scaleby;
    kx = 0.;
    sum += func(kx,ky,val);	//x DC term
    for (int j=1; j< N/2 ; j++) {
      kx = j*dk;
      val = *(zptr++) * scaleby;
      sum += func(kx,ky,val);
      sum += func(-kx,-ky,conj(val));
    }
    kx = dk*N/2;
    val = *(zptr++) * scaleby;
    sum += func(kx,ky,val);
  }

  // wrap to the negative ky's
  for (int i=-N/2; i< 0; i++) {
    ky = i*dk;
    val = *(zptr++) * scaleby;
    kx = 0.;
    sum += func(kx,ky,val);	//x DC term
    for (int j=1; j< N/2 ; j++) {
      kx = j*dk;
      val = *(zptr++) * scaleby;
      sum += func(kx,ky,val);
      sum += func(-kx,-ky,conj(val));
    }
    kx = dk*N/2;
    val = *(zptr++) * scaleby;
    sum += func(kx,ky,val);
  }
  sum *= dk*dk;
  return sum;
}

// Integrate kTable over d^2k (sum of all pixels * dk * dk)
DComplex
kTable::integratePixels() const {
  DComplex sum=0.;
  DComplex *zptr=array;
  // Do the positive y frequencies
  for (int i=0; i<= N/2; i++) {
    sum += *(zptr++) * scaleby;    // x DC term
    for (int j=1; j< N/2 ; j++) {
      sum += *(zptr) * scaleby;
      sum += conj(*(zptr++) * scaleby);
    }
    sum += *(zptr++) * scaleby;
  }
  // wrap to the negative ky's
  for (int i=-N/2; i< 0; i++) {
    sum += *(zptr++) * scaleby;    // x DC term
    for (int j=1; j< N/2 ; j++) {
      sum += *(zptr) * scaleby;
      sum += conj(*(zptr++) * scaleby);
    }
    sum += *(zptr++) * scaleby;
  }
  sum *= dk*dk;
  return sum;
}

// Make a new table that is function of old.
kTable*
kTable::Function( DComplex func(const double kx, 
				const double ky, 
				const DComplex val)) const {
  kTable *lhs = new kTable(N,dk);
  DComplex val;
  double kx, ky;
  DComplex *zptr=array;
  DComplex *lptr=lhs->array;
  // Do the positive y frequencies
  for (int i=0; i< N/2; i++) {
    ky = i*dk;
    for (int j=0; j<= N/2 ; j++) {
      kx = j*dk;
      val = *(zptr++) * scaleby;
      *(lptr++)= func(kx,ky,val);
    }
  }
  // wrap to the negative ky's
  for (int i=-N/2; i< 0; i++) {
    ky = i*dk;
    for (int j=0; j<= N/2 ; j++) {
      kx = j*dk;
      val = *(zptr++) * scaleby;
      *(lptr++)= func(kx,ky,val);
    }
  }
  return lhs;
}

// Transform to a single x point:
double
kTable::xval(double x, double y) const {
  // ??? check this:  don't evaluate if x not in fundamental period:
  x*=dk; y*=dk;
  if (x > 2*PI || y > 2*PI) throw FFTOutofRange();
  DComplex I(0.,1.);
  DComplex dxphase=exp(I*x);
  DComplex dyphase=exp(I*y);
  DComplex phase(1.,0.);
  DComplex z;
  double sum=0.;
  // y DC terms first:
  DComplex *zptr=array;
  // Do the positive y frequencies
  DComplex yphase=1.;
  for (int i=0; i< N/2; i++) {
    phase = yphase;
    z= *(zptr++);
    sum += (phase*z).real();	//x DC term
    for (int j=1; j< N/2 ; j++) {
      phase *= dxphase;
      z= *(zptr++);
      sum += (phase*z).real() * 2.;
    }
    phase *= dxphase;		//j=N/2 has no mirror:
    z= *(zptr++);
    sum += (phase*z).real();
    yphase *= dyphase;
  }

  // wrap to the negative ky's
  yphase = exp(I*(y*(-N/2)));
  for (int i=-N/2; i< 0; i++) {
    phase = yphase;
    z= *(zptr++);
    sum += (phase*z).real() * 2.;
    for (int j=1; j< N/2 ; j++) {
      phase *= dxphase;
      z= *(zptr++);
      sum += (phase*z).real() * 2.;
    }
    phase *= dxphase;		//j=N/2 has no mirror:
    z= *(zptr++);
    sum += (phase*z).real();
    yphase *= dyphase;
  }

  sum *= dk*dk*scaleby/(4.*PI*PI);	//inverse xform has 2pi in it.
  return sum;
}

void
// Translate the PSF to be for source at (x0,y0);
// ??need a sign flip here?
kTable::Translate(double x0, double y0) {
  // convert to phases:
  x0*=dk; y0*=dk;
  // too big will just be wrapping around:
  if (x0 > PI || y0 > PI) throw FFTOutofRange();
  DComplex I(0.,1.);
  DComplex dxphase=exp(DComplex(0.,x0));
  DComplex dyphase=exp(DComplex(0.,y0));
  DComplex phase(1.,0.);

  DComplex yphase=1.;
  DComplex z;

  DComplex *zptr=array;

  for (int i=0; i< N/2; i++) {
    phase = yphase;
    for (int j=0; j<= N/2 ; j++) {
      z = *zptr;
      *zptr = phase * z;
      phase *= dxphase;
      zptr++;
    }
    yphase *= dyphase;
  }

  // wrap to the negative ky's
  yphase = exp(I*((-N/2)*y0));
  for (int i=-N/2; i< 0; i++) {
    phase = yphase;
    for (int j=0; j<= N/2 ; j++) {
      z = *zptr;
      *zptr = phase* z;
      phase *= dxphase;
      zptr++;
    }
    yphase *= dyphase;
  }
  return;
}

      

//--------------------------------------------------
xTable::xTable(int _N, 
	       double _dx, 
	       double value) {
  N = 2*((_N+1)/2);	//Round size up to even.
  dx = _dx;
  get_array(value);
  valid = true;
  scaleby=1.;
  return;
}

size_t
xTable::index(int i, int j) const {
  // origin will be in center.
  i += N/2;
  j += N/2;
  if (i<0 || i>=N || j<0 || j>=N) throw FFTOutofRange() ;
  return i*N+j;
}

double xTable::xval(int i, int j) const {
  validCheck();
  return scaleby*array[index(i,j)];
}

void
xTable::xSet(int i, int j, double value) {
  validCheck();
  array[index(i,j)]=value/scaleby;
  return;
}

void
xTable::get_array(const double value) {
  array = (double*) fftw_malloc(sizeof(double)*N*N);
  for (int i=0; i<N*N; i++)
    array[i]=value;
  return;
}

void
xTable::copy_array(const xTable& rhs) {
  if (array!=0 && N!=rhs.N) {
    fftw_free(array);
    array = 0;
  }
  if (!rhs.valid || rhs.array==0) {
    valid=false;
    fftw_free(array);
    array = 0;
    return;
  } else
    valid = true;
  if (array==0)   array = (double*) fftw_malloc(sizeof(double)*N*N);
  for (int i=0; i<N*N; i++)
    array[i]=rhs.array[i];
  return;
}

void
xTable::kill_array() {
  fftw_free(array);
  array=0;
  return;
}

// Interpolate table (linearly) to some specific k:
double
xTable::xval(double x, double y) const {
  x /= dx;
  y /= dx;
  int i = static_cast<int> (floor(y));
  double fy = y-i;
  int j = static_cast<int> (floor(x));
  double fx = x-j;

  return (1-fx)*( (1-fy)*xval(i,j) + fy*xval(i+1,j)) +
    fx*( (1-fy)*xval(i,j+1) + fy*xval(i+1,j+1));
}

// Set table point nearest to some k
void
xTable::xSet(double x, double y, double value) {
  x /= dx;
  y /= dx;
  int j = static_cast<int> (floor(x+0.5));
  int i = static_cast<int> (floor(y+0.5));
  xSet(i,j,value);
  return;
}
  
// Fill table from a function:
void
xTable::fill( double func(const double x, const double y)) {
  double *zptr=array;
  double x, y;
  for (int i=0; i<N; i++) {
   y = (i-N/2)*dx;
   for (int j=0; j< N ; j++) {
     x = (j-N/2)*dx;
     *(zptr++) = func(x,y);
   }
  }
  return;
}

// Integrate a function over x - can be function of x or of PSF(x)
// Setting the Boolean flag gives sum over samples, not integral.
double
xTable::Integrate( double func(const double x, 
			       const double y, 
			       const double val),
		   bool  sumonly) const {
  double sum=0.;
  double val;
  double x, y;
  double *zptr=array;

  for (int i=0; i< N; i++) {
    y = (i-N/2)*dx;
    for (int j=0; j< N ; j++) {
      x = (j-N/2)*dx;
      val = *(zptr++) * scaleby;
      sum += func(x,y,val);
    }
  }

  if (!sumonly) sum *= dx*dx;
  return sum;
}

double
xTable::integratePixels() const {
  double sum=0.;
  double *zptr=array;
  for (int i=-N/2; i< N/2; i++) 
    for (int j=-N/2; j< N/2; j++) {
      sum += *(zptr++) * scaleby;
    }
  sum *= dx*dx;
  return (double) sum;
}

// Transform to a single k point:
DComplex
xTable::kval(double kx, double ky) const {
  // check this:  don't evaluate if x not in fundamental period:
  kx*=dx; ky*=dx;
  if (kx > 2*PI || ky > 2*PI) throw FFTOutofRange();
  DComplex I(0.,1.);
  DComplex dxphase=exp(-I*kx);
  DComplex dyphase=exp(-I*ky);
  DComplex phase(1.,0.);
  DComplex z;
  DComplex sum=0.;

  double *zptr=array;
  DComplex yphase=exp(I*(ky*N/2));
  for (int i=0; i< N; i++) {
    phase = yphase;
    phase *= exp(I*(kx*N/2));
    for (int j=0; j< N ; j++) {
      sum += phase* (*(zptr++));
      phase *= dxphase;
    }
    yphase *= dyphase;
  }
  sum *= dx*dx*scaleby;
  return sum;
}

// Have FFTW develop "wisdom" on doing this kind of transform
void 
kTable::fftwMeasure() const {
  DComplex* t_array = 
    (DComplex*) fftw_malloc(sizeof(DComplex)*N*(N/2+1));
  // Copy data into new array to avoid NaN's, etc., but not bothering
  // with scaling, etc.
  for (int i=0; i<N*(N/2+1); i++)
    t_array[i] = array[i];

  xTable *xt = new xTable( N, 2*PI/(N*dk) );

  fftw_plan plan = 
    fftw_plan_dft_c2r_2d(N, N,
			 reinterpret_cast<fftw_complex*> (t_array), xt->array,
			 FFTW_MEASURE);
  if (plan==NULL) throw FFTInvalid();
  delete xt;
  fftw_free(t_array);
  fftw_destroy_plan(plan);
}

// Fourier transform from (complex) k to x:
xTable*
kTable::Transform() const {
  // We'll need a new k array because FFTW kills the k array in this
  // operation.  Also, to put x=0 in center of array, we need to flop
  // every other sign of k array, and need to scale.

  DComplex* t_array = 
    (DComplex*) fftw_malloc(sizeof(DComplex)*N*(N/2+1));
  double fac = scaleby * dk * dk / (4*PI*PI);
  long int ind=0;
  for (int i=0; i<N; i++)
    for (int j=0; j<=N/2; j++) {
      if ( (i+j)%2==0) t_array[ind]=fac * array[ind];
      else t_array[ind] = -fac* array[ind];
      ind++;
    }

  xTable *xt = new xTable( N, 2*PI/(N*dk) );

  fftw_plan plan = 
    fftw_plan_dft_c2r_2d(N, N,
			 reinterpret_cast<fftw_complex*> (t_array), xt->array,
			 FFTW_ESTIMATE);
  if (plan==NULL) throw FFTInvalid();

  // Run the transform:
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  fftw_free(t_array);
  return xt;
}

// same function, takes xTable reference as agrument 
void
kTable::Transform(xTable& xt) const {

  // ??? check proper dimensions for xt ?
  Assert(N==xt.getN());

  // We'll need a new k array because FFTW kills the k array in this
  // operation.  Also, to put x=0 in center of array, we need to flop
  // every other sign of k array, and need to scale.

  DComplex* t_array = 
    (DComplex*) fftw_malloc(sizeof(DComplex)*N*(N/2+1));
  double fac = scaleby * dk * dk / (4*PI*PI);
  long int ind=0;
  for (int i=0; i<N; i++)
    for (int j=0; j<=N/2; j++) {
      if ( (i+j)%2==0) t_array[ind]=fac * array[ind];
      else t_array[ind] = -fac* array[ind];
      ind++;
    }

  fftw_plan plan = 
    fftw_plan_dft_c2r_2d(N, N,
			 reinterpret_cast<fftw_complex*> (t_array), xt.array,
			 FFTW_ESTIMATE);
  if (plan==NULL) throw FFTInvalid();

  // Run the transform:
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  fftw_free(t_array);

  xt.dx = 2*PI/(N*dk);
  xt.scaleby = 1.;
  xt.valid = true;
}

void 
xTable::fftwMeasure() const {
  // Make a new copy of data array since measurement will overwrite:
  double* t_array = 
    (double*) fftw_malloc(sizeof(double)*N*N);
  // Copy data into new array to avoid NaN's, etc., but not bothering
  // with scaling, etc.
  for (int i=0; i<N*N; i++)
    t_array[i] = array[i];

  kTable *kt = new kTable( N, 2*PI/(N*dx) );

  fftw_plan plan = 
    fftw_plan_dft_r2c_2d(N,N,
			 t_array, reinterpret_cast<fftw_complex*> (kt->array),
			 FFTW_MEASURE);
  if (plan==NULL) throw FFTInvalid();

  delete kt;
  fftw_free(t_array);
  fftw_destroy_plan(plan);
}

// Fourier transform from x back to (complex) k:
kTable*
xTable::Transform() const {

  kTable *kt = new kTable( N, 2*PI/(N*dx) );

  fftw_plan plan = 
    fftw_plan_dft_r2c_2d(N,N,
			 array, reinterpret_cast<fftw_complex*> (kt->array),
			 FFTW_ESTIMATE);
  if (plan==NULL) throw FFTInvalid();
  fftw_execute(plan);
  fftw_destroy_plan(plan);

  // Now scale the k spectrum and flip signs for x=0 in middle.
  double fac = scaleby * dx * dx; 
  size_t ind=0;
  for (int i=0; i<N; i++)
    for (int j=0; j<=N/2; j++) {
      if ( (i+j)%2==0) kt->array[ind] *= fac;
      else kt->array[ind] *= -fac;
      ind++;
    }

  return kt;
}
