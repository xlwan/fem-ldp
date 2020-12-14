/**************************************************************************
 This class provides utilities of 1D smoother through projection or 
 elliptic operator. 

 Author: Xiaoliang Wan
 Last updated: 06/02/2014
 *************************************************************************/

#ifndef SMOOTHER1D_H_
#define SMOOTHER1D_H_

class Smoother1D{
 public:
  int     nTE; // number of elements.

  double*  dT; // element size.

  double*  xb; // node basis.
  double*  rb; // rhs for xb.

  // matrices for boundary modes
  double**  Q; // tridiagonal mass matrix.
  double**  A; // tridiagonal stiffness matrix for xb.

  Smoother1D();
  virtual ~Smoother1D();
  
  void form_smoother_operator(int nTE_in, double* mesh);

  // smooth piecewise constant
  void smooth_p0_to_p1_projection(double* p0, double* p1);

  // smooth piecewise linear element
  void smooth_p1_to_p1_elliptic_cg(double* p1_before, double* p1_after);
  void smooth_p1_to_p1_elliptic_jacobi(double* p1_before, double* p1_after);
  void smooth_p1_to_p1_elliptic_gs(double* p1_before, double* p1_after);

 private:
  void Axv(double* x_in, double* x_out);
  void Rxv(double* x_in, double* x_out); // R is A without diagonal elements.
};

#endif
