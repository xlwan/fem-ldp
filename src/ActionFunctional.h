/**************************************************************************
 * ActionFunctional.h : Class ActionFunctional
 *
 * ActionFunctional provides a general definition of the approximation 
 * of the path, and an interface between the nonlinear conjugate gradient 
 * solver and the elementwise computation of action functional and its 
 * gradient for a specific problem. 
 *
 * Author: Xiaoliang Wan
 **************************************************************************/

#ifndef ACTIONFUNCTIONAL_H_
#define ACTIONFUNCTIONAL_H_

#include "CG_Descent.h"
#include "ElementT.h"
#include "HPadaptivity.h"

class ActionFunctional //: public CG_Descent 
{
 public:
  int       dim;   // dimension
  int       nTE;   // number of elements in time direction.

  double      T;   // integration time: [0,T].
  int        nG;   // total number of unknowns.

  // Boundary conditions of MAP at t= 0 & T.
  double*    J0;   // expansion coefficients of the end of the path at t=0.
  double*    JT;   // expansion coefficients of the end of the path at t=T.

  // Rewritten MAP as a list of finite elements. 
  double*    Tmesh;  // temporal discretization.  
  int*       nT;     // number of modes in each element.

  ElementT*  mapList_head; // head of the element list for map.
  ElementT** mapTable;     // organize the elements in a table. 

  // need something for adaptivity.
  int nTE_adpt;  // element number for the new mesh. 

  Smoother1D* gRe; // the 1D smoother for a posteriori error estimate.

  ElementT** mapTable_adpt;     // organize the elements in a table.

  ElementT** sort_elmt_list; // sort all elements;
  // the dimensionality of these two arrays will be determined according to 
  // the current status of the partition.

  // For hp adaptivity
  HP_adaptivity* hpAdpt;

  // optimization solver.
  CG_Descent* cgDescent;

  ActionFunctional(int dim_in);
  virtual ~ActionFunctional();

  /************ main functions for the adaptive MAM *****************/
  // nonlinear CG descent solver.
  void prepare_CG_Solver(cg_parameter* cgparam);
  void CG_Descent_Solver();

  // compute the action functional.
  double compute_AF  (double* path); 

  // compute the gradient of the action functional.
  void   compute_Grad (double* path, double* grad); 

  // compute both action functional and gradient.
  double compute_Grad_AF(double* path, double* grad);
  /********************************************************************/


  /************ generation and initialization of path *****************/
  // set B.C.s of the MAP.
  void setBCs(double* J0_in, double* JT_in);

  // generate element list.
  void generate_Elmt_List(FILE* fp);
  void generate_Elmt_List(int nT_in, int nTE_in);

  // generate global map.
  void form_gmap();

  // set path from a given function.
  void func2path();
  /********************************************************************/

  /************* switch between a vector and element representation of the path ***********/
  // scatter the MAP to each element.
  void scatter_Path(double* path);

  // gather the path from each element.
  void gather_Path(double* path);

  // gather the gradient of the action functional from each element.
  void gather_Grad(double* path);
  /****************************************************************************************/

  /********** hp adaptivity utilities  ******************
   *
   * All specific routines will be found in the class object hpAdpt.
   *
   ******************************************************/
  // generate the adaptive mesh
  void generate_Adaptive_Mesh();

  // make appropriate adjustment for the new mesh
  void post_adpt_adjustment();

  /*************** wrappers for cg solver *******************************************
   * The CG solver requires the function pointers to the functions of 
   * computing the action functional, gradient or both. 
   *
   * The argument void* instance needs to be trasfered to an ActionFunctional pointer.
   **********************************************************************************/
  static double Get_Action   (void* instance, double* x, int n);
  static void   Get_Gradient (void* instance, double* g, double* x, int n);
  static double Get_Act_Grad (void* instance, double* g, double* x, int n);

  /********************* io routines ***********************/
  // backup the MAP very a certain number of iteration steps.
  int  backup_MAP(const char* path1);

  // dump the information of MAP and time mesh
  void dump_MAP(double* x);

  double summary_of_initial_path(double* x, double* g);

  void add_info_to_report(FILE* fp);
  /**********************************************************/

 private:
  double f_start;       // initial value of the action functional.
  double gnorm_start;   // initial Linfy norm of the gradient.
  double gnorm2_start;  // initial L2 norm of the gradient.

  Param* param_restart;

  int ifixT;      // integration time is fixed or not.
};

#endif /* ACTIONFUNCTIONAL_H_ */
