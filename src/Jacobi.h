/*
 * Jacobi.h
 * Basic utilities related to Jacobi polynomials on the interval [-1,1].
 *
 * Author: Xiaoliang Wan
 */

#ifndef JACOBI_H_
#define JACOBI_H_

typedef enum gType{
  Gauss,    // Gauss quadrature points.
  Lobatto,  // Gauss-Lobatto quadrature points.
  RadauL,   // Gauss-Radau quadrture points including the left end.
  RadauR    // Gauss-Radau quadrture points including the right end.
}gType;

typedef struct qstore{
  double  alpha;      /* (1-x)^alpha*(1+x)^beta */
  double  beta;       /* (1-x)^alpha*(1+x)^beta */
  int     n;          /* quadrature order    */
  double* qpt;        /* quadrature points   */
  double* wgt;        /* integration weights */
  gType   type;       /* type of quadrature points */
  struct qstore *next;
} Qstore;

class Jacobi {
  double alpha;
  double beta;
 public:
  static Qstore* qw;
  
  // get the Gauss-type quadrature points (qpt) and the corresponding weights (wgt).
  void getQW(int npt, double** qpt, double** wgt, gType type);
  // npt: number of points.
  // qpt: quadrature points in [-1,1].
  // wgt: integration weights.
  
  // compute the values (val) and first-order derivatives (div) of n-th order Jacobi polynomrilas at x.
  void getVD(int npt, double* x, double* val, double* div, int n);
  // npt: number of points.
  //   x: arbitrarily given points in [-1,1].
  // val: values of Jacobi polynomial at x.
  // div: first-order derivative of Jacobi polynomial at x.
  //   n: order of the Jacobi polynomial.
  
  // define the derivative matrix D and its transpose Dt on npt Gauss-type quadrature points qpt.
  void getDmatrix(double* D, double* Dt, int npt, double* qpt, gType type);
  //   D: derivative matrix D.
  //  Dt: transpose of D.
  // npt: number of quadrature points.
  // qpt: quadrature points.
  //type: type of quadrature points.
  
  // define the interpolation operators through Lagrangian interpolants from npt quadrature points to
  // nx arbitrarily given points x.
  void getImatrix(double *Q2X, int npt, double *qpt, int nx, double *x, gType type);
  // Q2X: interpolation matrix.
  // npt: number of points.
  // qpt: quadrature points.
  //  nx: number of points x.
  //   x: arbitrarily given points x.
  //type: type of Gauss quadrature points.
  
  Jacobi(double alpha_in, double beta_in);
  virtual ~Jacobi();
};


#endif /* JACOBI_H_ */
