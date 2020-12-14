/*
 * Basis.h
 *
 * class basis provides the information for one-dimensional basis.
 *  
 * class Basis provides the two-dimensional basis for plane parallel flow. 
 * Author: Xiaoliang Wan
 */

#ifndef BASIS_H_
#define BASIS_H_

#include "Jacobi.h"
#include "fftw3.h"

/*
  basis defines a one-dimensional approximation basis. 
  It includes three types:
  'v': P_n-P_{n+2} with P_n being Legendre polynomials. 
       All basis functions are equal to zero at the boundary -1 and 1.
  'p': P_n. The regular Legendre polynomials.
  't': Spectral/hp element basis. See Karniadakis and Sherwin's book.
 */

class basis{
 public:
  char     bt;   // a certain type of basis, e.g. a certain combination of Legendre polynomials.
  int      nbm;  // number of basis modes.
  int      npt;  // number of quadrature points.
  gType    gt;   // type of quadrature points.
  double*  val;  // function values on the quadrature points.
  double*  D1;   // first-order derivative on quadrature points.
  double*  D2;   // second-order derivative on quadrature points.
  double*  Dm;   // differentiation matrix.
  double*  DmT;  // transpose of the differentiation matrix.
  int      nband;// if the mass matrix is a symmetric banded one.
  double*  mass; // mass matrix.
  double*  invm; // Cholesky factoraization of the mass matrix.
  double*  wint;  // integration weights.

  double*  coef_orth; // projection to previous modes

  double* work;

  basis*   next;

  basis(char bt_in, int nbm_in, int npt_in, gType gt_in);
  virtual ~basis();

  // set up the basis.
  void setbasis();

  // compute the value and derivative of a certain basis mode.
  void getVD(int n, double* x, double* p, double* dp, int mode);

  // transform between modal and physical spaces.
  void J2Q(double* J, double* Q);
  void Q2J(double* Q, double* J);

  // gradient
  void Grad_J2Q(double* J,  double* Q);
  void Grad_Q2Q(double* Qf, double* Qg);

  // projection of the mode of highest order to previous modes.
  void getOrthProjection(int m, double* coef);
  void setOrthProjection();

  // (nbm-1)th derivative of the mode of highest order.
  void getD2Cnst(double* v, double elmt_size);
  double getD2Linear(double elmt_size);
};

class Basis { 
public:
  int     nX;      // number of Fourier modes in x direction.
  int     nY;      // number of basis modes in y direction. (P_n-P_{n+2})
	
  // Information for the x direction: Fourier transform is assumed for the x direction.
  int     nptX;
  double* bvX;             // basis values on quadrature points.
  double* d1X;             // first-order derivative on quadrature points.
  double* d2X;             // second-order derivative on quadrature points.

  // Information for the y direction: Legendre polynomials or their combinations.
  int     nptY;            // number of quadrature points.
  char    btypeY;          // basis type. could be different for velocity and pressure fields.
  gType   gtypeY;          // type of quadrature points.
  double* bvY;             // basis values on quadrature points.
  double* d1Y;             // first-order derivative on quadrature points.
  double* d2Y;             // second-order derivative on quadrature points.
  double* DY;              // differential matrix defined by the basis in y direction.
  double* DtY;             // transpose of the differentiation matrix.

  double* wgD;             // integration weights.

  // Information for the mass matrix.
  int     kd;              // sub- or superdiagonals of the mass matrix.
  double* massY_Cholesky;  // Cholesky factors of mass matrix of basis in the y direction.

  // Information for FFTW.
  fftw_r2r_kind kindF; // Forward FFT sign:  R2HC
  fftw_r2r_kind kindB; // Backward FFT sign: HC2R
  fftw_plan     pF2B;  // Forward  FFT plan
  fftw_plan     pB2F;  // Backward FFT plan

  // working space
  double *wkF;
  double *wkB;
  double *wkD;
  
  Basis(int nX_in, int nY_in, int nptX_in, int nptY_in, char btypeY_in, gType gtypeY_in);
  virtual ~Basis();
  
  void setBasis(); // initialize the basis.
  
  void J2Q(double* J, double* Q); // modal space to physical space.
  void Q2J(double* Q, double* J); // physical space to modal space.
  
  void Grad_J2Q(double* J, double* dxQ, double* dyQ); 
  // Gradient, starting from modal space, end up with physical space.
  
  void Grad_Q2Q(double* Q, double* dxQ, double* dyQ); 
  // Gradient, starting from physical space, end up with physical space.
};

#endif /* BASIS_H_ */
