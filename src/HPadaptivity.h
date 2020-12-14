/**************************************************************************
 * HP_adaptivity provides a posteriori error estimate of the approximated 
 * transition path using the derivative recovery technique. 
 *
 * Smoother1D provides utilities of 1D smoother through projection or 
 * elliptic operator. 
 *
 * Author: Xiaoliang Wan
 *************************************************************************/

#ifndef HP_ADAPTIVITY_H_
#define HP_ADAPTIVITY_H_

#include "ElementT.h"

class Smoother1D{
 public:
  int      nE; // number of elements.

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
  void smooth_p1_to_p1_elliptic_cg     (double* p1_before, double* p1_after);
  void smooth_p1_to_p1_elliptic_jacobi (double* p1_before, double* p1_after);
  void smooth_p1_to_p1_elliptic_gs     (double* p1_before, double* p1_after);

 private:
  void Axv(double* x_in, double* x_out);
  void Rxv(double* x_in, double* x_out); // R is A without diagonal elements.
};

class HP_adaptivity{
 public:
  int                dim;  // number of dimension.
  
  int          nE_before;  // number of elements before adaptivity.
  int        nDOF_before;  // number of DOFs before adaptivity
  ElementT** elmt_before;  // element list before adaptivity.

  int           nE_after;  // number of elements after adaptivity
  int         nDOF_after;  // number of DOFs after adaptivity
  ElementT**  elmt_after;  // element list after adaptivity.
  
  int         nDOF_model;  // number of DOFs for model approximation

  Smoother1D* devRecover; // derivative recovery technique.

  HP_adaptivity(int dim_in, int nE_in, ElementT** elmt_old_in);
  virtual ~HP_adaptivity();

  // prepare the adaptivity.
  void set_param();
  void pre_Adpt();
  void a_posteriori_error_indicators();
  void compute_theta(); // compute the error indicator theta.
  void compute_eta();   // compute the local error indicator.

  // deal with model approximation using h-refinement
  void update_Adpt_ElmtList_Model();
  int  necessity_for_model_adpt();


  // deal with numerical approximation using hp-refinement.
  void update_Adpt_ElmtList_Approximation();
  ElementT* hp_refinement(ElementT* EL);


  // reshuffle the element list according to eta.
  void reshuffle_eta  (ElementT** a, int N, ElementT* EL);
  void reshuffle_eta  (ElementT** a, int N, ElementT* EL, ElementT* ER);
  void reshuffle_theta(ElementT** a, int N, ElementT* EL, ElementT* ER);
  
 private:
  int pmax; // highest poly order
  int pmin; // lowest poly order

  // sort the elements.
  void swap_E(ElementT** left, ElementT** right);
  int  partition_theta_descend(ElementT** a, int left, int right);
  int  partition_eta_descend(ElementT** a, int left, int right);
  void quicksort_theta_descend(ElementT** a, int left, int right);
  void quicksort_eta_descend(ElementT** a, int left, int right);
  int  iwhere_eta_descend(ElementT** a, int lb, int hb, double target);
  int  iwhere_theta_descend(ElementT** a, int lb, int hb, double target);

  double ave;      // average alpha.
};

#endif
