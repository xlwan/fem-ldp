/*********************************************************************
 * ElementT.h: class ElementT
 *
 * The class ElementT provides a finite element representation of the 
 * minimal action path. Most of the operations, such as computation of
 * action functional and  its gradient, will be decomposed into local 
 * operations in each element.
 * 
 * Author: Xiaoliang Wan
 *********************************************************************/

#ifndef ELEMENTT_H_
#define ELEMENTT_H_

#include "Basis.h"

class ElementT{
 public:  
  /*** basic approximation information *******************/
  double    dT; // time step size.
  double    lb; // lower bound. [lb, hb]
  double    hb; // higher bound. [lb, hb]
  int       nT; // number of basis modes in T direction.
  int     nptT; // number of quadrature points for T.
  /*******************************************************/

  double     T; // integration time.

  /*** the following several variables will be determined by a certain problem ***/
  int      dim; // dimension of Ju.
  double**  Ju; // coefficients of solution u.
  double**  Qu; // values on quadrature points;
  double**  Du; // derivative. 

  double**  fut; // force term.
  double**  fut_d; // dx
  double**  fut_b; // b(x)
  double*** fgt; // grad of force term.

  int* loc_mapping; 
  /*******************************************************************************/

  basis* bt; // basis for approximation in time.

  // we define a bio-directional list.
  ElementT* next;
  

  // indicators for adaptivity.
  double   theta;         // local indicator for arc length constraint.
  double   eta;           // local error indicator.
  double   alpha;         // local regularity indicator. 
  double** p_extension;  // coeficient of the (p+1)th order term.
  double** p_ext_GLQ;     // (p+1)th order extension at Gauss Lobatto quadrature points.

  // this is a list of approximation basis on a reference element.
  static basis* blist; // list of finite element basis.

  // look for information of a certain type of approximation basis.
  static basis* getbasis(int nx, int nptx, char uvp);


  ElementT(int nT_in, int dim_in);
  virtual ~ElementT();

  /***********  preparation of the element ***************/
  // geometry info
  void set_geo_info(double lb_in, double hb_in);

  // allocate the memory.
 // this function allocate the memory for the path.
  void mem_alloc_map();  
  /*******************************************************/


  /****** local contribution to the action functional and gradient *******/  
  // prepare the computation of action functional.
  void pre_compt_AF();

  // local contribution  for the evalution of the action functional.
  double local_val_AF();

  // prepare the computation of gradient.
  void pre_compt_Grad();
 
  // local contribution for the gradient of the action functional.
  void local_grad_AF();
  double delta_grad(int iT, int jU);

  // local contribution of <dx, dx>
  double local_vd();
  
  // local contribution of <b(x), b(x)>
  double local_vb();

  // local arc length constraint
  double local_arc_length_constraint();
  /************************************************************************/

  /*********** following are utilities for adaptivity *********************/
  double get_arc_length(double t);

  // compute the pth order derivative of approximate solution.
  void compt_derivative_order_p(int p, double* vc);

  // h_refinement
  ElementT* h_refinement();

  // p_refinement
  void p_refinement();

  // compute error indicator.
  void compute_eta();

  // compute theta indicator
  void compute_theta();

  // relocate memory for p refinement.
  void relocate_mem_for_p_refinement();
  /************************************************************************/

  /************ following are some useful utilities needed *****************/
  // compute the path on quadrature points.
  void J2Q();
  void Q2J();

  // snapshot of the path for a certain quadrautre point or time.
  void snapshot(int in, double* p, double* pt);
  void snapshot(double t, double* v, double* vt);

  // initialize the path Ju from some given functions.
  void func2Ju(double (**f)(double));

 private:
  double* Imat_L; // interpolation matrix for the left element;
  double* Imat_R; // interpolation matrix for the right element;
};

#endif /* ELEMENTT_H */
