/***************************
 * Author: Xiaoliang Wan
 ***************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include "FortranMapping.h"
#include "Common.h"
#include "HPadaptivity.h"

Smoother1D::Smoother1D()
{
 
  dT = 0;

  xb = 0;
  rb = 0;
  Q  = 0;
  A  = 0;
}

Smoother1D::~Smoother1D()
{
  if(dT) delete[] dT;
  if(xb) delete[] xb;
  if(rb) delete[] rb;
  if(Q)  delete_matrix(Q);
  if(A)  delete_matrix(A);
}

void Smoother1D::form_smoother_operator(int nE_in, double* mesh)
{
  nE = nE_in;


  int n = nE+1;
  if(Q)
    {
      delete_matrix(Q);
      Q = dmatrix(3, n);
    }
  else
    Q = dmatrix(3, n);

  if(A)
    {
      delete_matrix(A);
      A = dmatrix(3, n);
    }
  else
    A = dmatrix(3, n);

  if(xb)
    {
      delete[] xb;
      xb = new double[n];
    }
  else
    xb = new double[n];

  if(rb)
    {
      delete[] rb;
      rb = new double[n];
    }
  else
    rb = new double[n];

  if(dT)
    {
      delete[] dT;
      dT = new double[nE];
    }
  else
    dT = new double[nE];

  for(int i = 0; i < nE; i++)
    dT[i] = mesh[i+1] - mesh[i];

  dzero(3*n, Q[0], 1);
  dzero(3*n, A[0], 1);

  // mass matrix - tridiagonal.
  for(int i = 1; i < n; i++)
    Q[0][i] = dT[i-1]/6.0;

  Q[1][0] = dT[0]/3.0;
  for(int i = 1; i < n-1; i++)
    Q[1][i] = (dT[i-1]+dT[i])/3.0;
  Q[1][n-1] = dT[nE-1]/3.0;

  for(int i = 0; i < (n-1); i++)
    Q[2][i] = dT[i]/6.0;
  
  // stiffness matrix - tridiagonal.
  for(int i = 1; i < n; i++)
    A[0][i] = -1.0/dT[i-1] + Q[0][i];

  A[1][0] = 1.0/dT[0] + Q[1][0];
  for(int i = 1; i < n-1; i++)
    A[1][i] = (1.0/dT[i-1] + 1.0/dT[i]) + Q[1][i];
  A[1][n-1] = 1.0/dT[nE-1] + Q[1][n-1];

  for(int i = 0; i < n-1; i++)
    A[2][i] = -1.0/dT[i] + Q[2][i];

  return;
}

//smooth piecewise constant - through projection
void Smoother1D::smooth_p0_to_p1_projection(double* p0, double* p1)
{
  int n = nE+1;

  rb[0] = p0[0]*dT[0]/2.0;
  for(int i = 1; i < n-1; i++)
    rb[i] = (p0[i-1]*dT[i-1] + p0[i]*dT[i])/2.0;
  rb[n-1] = p0[nE-1]*dT[nE-1]/2.0;

  tridag_solver(Q[0], Q[1], Q[2], rb, p1, n);

  return;
}

//smooth piecewise linear element using conjugate gradient
void Smoother1D::smooth_p1_to_p1_elliptic_cg(double* p1_before, double* p1_after)
{
  int n = nE+1;

  double* work = new double[4*n];

  double* p_next  = work;
  double* p_now   = p_next + n;
  double* ap_now  = p_now  + n;
  double* r_next  = ap_now + n;
  double* r_now   = rb;

  double* x_now   = xb;
  double* x_next  = p1_after;

  // x0
  dcopyC(n, p1_before, 1, x_now, 1);

  // r0
  Axv(x_now, r_now);
  dscalC(n, -1.0, r_now, 1);

  double tmp = sqrt(ddotC(n, r_now, 1, r_now, 1));
  if(tmp < 1e-12)
    {
      dcopyC(n, p1_before, 1, x_next, 1);
      return;
    }

  // p0
  dcopyC(n, r_now, 1, p_now, 1);

  double alpha;
  double beta;
  for(int i = 0; i < 1; i++)
    {
      // update x, r, p from _now to _next
      Axv(p_now, ap_now);
      alpha = ddotC(n, r_now, 1, r_now, 1)/ddotC(n, p_now, 1, ap_now, 1);

      dcopyC(n, x_now, 1, x_next, 1);
      daxpyC(n, alpha, p_now, 1, x_next, 1);

      dcopyC(n, r_now, 1, r_next, 1);
      daxpyC(n, -alpha, ap_now, 1, r_next, 1);

      beta = ddotC(n, r_next, 1, r_next, 1)/ddotC(n, r_now, 1, r_now, 1);
      
      dcopyC(n, r_next, 1, p_next, 1);
      daxpyC(n, beta, p_now, 1, p_next, 1);

      // prepare next iteration.
      dcopyC(n, x_next, 1, x_now, 1);
      Axv(x_now, r_now);
      dscalC(n, -1.0, r_now , 1);
      dcopyC(n, p_next, 1, p_now, 1);
    }
 
  delete[] work;
  return;
}

// smooth piecewise linear element using jacobi iteration
void Smoother1D::smooth_p1_to_p1_elliptic_jacobi(double* p1_before, double* p1_after)
{
  int n = nE+1;

  double* x_now  = xb;
  double* x_next = p1_after;

  double* work = rb;

  dcopyC(n, p1_before, 1, x_now, 1);

  for(int i = 0; i < 1; i++)
    {
      Rxv(x_now, work);
      dscalC(n, -1.0, work, 1);

      for(int j = 0; j < n; j++)
	x_next[j] = work[j]/A[1][j];

      dcopyC(n, x_next, 1, x_now, 1);
    }

  return;
}

// smooth piecewise liear element using Gauss-Seidel iteration
void Smoother1D::smooth_p1_to_p1_elliptic_gs(double* p1_before, double* p1_after)
{
  int n = nE+1;

  double* x_now  = xb;
  double* x_next = p1_after;

  dcopyC(n, p1_before, 1, x_now, 1);

  for(int i = 0; i < 1; i++)
    {
      x_now[0] = -(A[0][1]*x_now[1])/A[1][0];
      for(int j = 1; j < n-1; j++)
	x_now[j] = -(A[2][j-1]*x_now[j-1]+A[0][j+1]*x_now[j+1])/A[1][j];
      x_now[n-1]= -(A[2][n-2]*x_now[n-2])/A[1][n-1];

      dcopyC(n, x_now, 1, x_next, 1);
    }

  return;
}

// matrix vector multiplication: Ax.
void Smoother1D::Axv(double* x_in, double* x_out)
{
  int n = nE+1;

  x_out[0] = A[1][0]*x_in[0] + A[0][1]*x_in[1];
  for(int i = 1; i < n-1; i++)
    x_out[i] = A[2][i-1]*x_in[i-1]+A[1][i]*x_in[i]+A[0][i+1]*x_in[i+1];
  x_out[n-1]= A[2][n-2]*x_in[n-2] + A[1][n-1]*x_in[n-1];

  return;
}

// matrix vector multiplication: Rx - R is A without the diagonal line.
void Smoother1D::Rxv(double* x_in, double* x_out)
{
  int n = nE+1;

  x_out[0] = A[0][1]*x_in[1];
  for(int i = 1; i < n-1; i++)
    x_out[i] = A[2][i-1]*x_in[i-1]+A[0][i+1]*x_in[i+1];
  x_out[n-1]= A[2][n-2]*x_in[n-2];

  return;
}

/******************************************************************************/

/******************      subroutines of HP_adaptivity    **********************/

/******************************************************************************/

HP_adaptivity::HP_adaptivity(int dim_in, int nE_in, ElementT** elmt_list_in)
{
  dim = dim_in;

  nE_before   = nE_in;
  elmt_before = elmt_list_in;

#ifdef MPISRC
  nDOF_before = 0;
  for(int i = 0; i  < pllinfo.nprocs; i++)
    nDOF_before += pllinfo.nd[i];
#else
  nDOF_before = 0;
  for(int i = 0; i < (nE_before-1); i++)
    nDOF_before += (elmt_before[i]->nT - 1)*dim;
  nDOF_before += (elmt_before[nE_before-1]->nT - 2)*dim;
#endif
  
  nDOF_model = 0;

  devRecover = new Smoother1D();

  elmt_after = 0;
}

HP_adaptivity::~HP_adaptivity()
{
  if(elmt_after) delete[] elmt_after;
  if(devRecover) //devRecover->~Smoother1D();
    delete devRecover;
}

void HP_adaptivity::pre_Adpt()
{
  // prepare indicators.
  a_posteriori_error_indicators();

  // adaptivity starts from the current DOF and number of elements.
  nDOF_after = nDOF_before;
  nE_after   = nE_before;

  // compute the average of alpha
  ave = 0;
  for(int i = 0; i < nE_before; i++)
    ave += elmt_before[i]->alpha;

#ifdef MPISRC
  double tmp = 0;
  gdsum(&ave, 1, &tmp);

  ave /= (double)pllinfo.gne;
#else
  ave /= (double)nE_before;
#endif

  // resize elmt_after.
  // The doubled number of elements is more than enough for the adaptive mesh.
  if(elmt_after) delete[] elmt_after;
  elmt_after = new ElementT*[3*nE_before];

  for(int i = 0; i < 3*nE_before; i++) 
    elmt_after[i] = 0;

  // the adaptivity starts from theta.
  for(int i = 0; i < nE_before; i++)
    elmt_after[i] = elmt_before[i];

  printf("ave is %20.14e with nTE %d\n", ave, nE_before);
  
  return;
}


// deal with numerical approximation
void HP_adaptivity::update_Adpt_ElmtList_Approximation()
{
  ElementT* EL;
  ElementT* ER;
 
  // sort eta.
  quicksort_eta_descend(elmt_after, 0, nE_after-1);

  int istrategy = (int)get_param("is_strategy_one");
  
  int itmp;
  double tmp;
  double eta_max;

  // istrategy = 1 means that we do adaptivity according to DOF.
  if(istrategy == 1)
    {  
      // compute the average polynomial order.
      double pve = 0;
      for(int i = 0; i < nE_before; i++)
	pve += (elmt_before[i]->nT - 1);
      
#ifdef MPISRC
      tmp = 0;
      gdsum(&pve, 1, &tmp);
      pve = pve/(double)pllinfo.gne;
#else
      pve = pve/(double)nE_before;
#endif     
 
      // compute the desired number of degrees of freedom.
      int nDOF_target = (int) ceil ((double)nDOF_before * (1 + pve) / pve);
      
      printf("pve is %lf, nDOF_before is %d, nDOF_target is %d\n", pve, nDOF_before, nDOF_target);

      while(nDOF_after < nDOF_target)
	{
#ifdef MPISRC
          tmp = DBL_MIN;
          eta_max = elmt_after[0]->eta;
          gdmax(&eta_max, 1, &tmp);

          //printf("eta_max: %20.14e; first: %20.14e\n", eta_max, elmt_after[0]->eta);
          //exit_code(-1);

          if(eta_max == elmt_after[0]->eta)
            {
              EL = elmt_after[0];
              ER = hp_refinement(EL);

              if(!ER)
                {
                  // p refinement
                  reshuffle_eta(elmt_after, nE_after, EL);
                  nDOF_after += dim;
                }
              else
                {
                  // h refinement
                  reshuffle_eta(elmt_after, nE_after, EL, ER);
                  nDOF_after += (EL->nT - 1)*dim;
                  nE_after++;
                }
            }
          
          itmp = INT_MIN;
          gimax(&nDOF_after, 1, &itmp);

          //ROOTONLY
          //  printf("current DOFs: %d; target DOFs: %d\n", nDOF_after, nDOF_target);
#else
	  EL = elmt_after[0];	  
	  ER = hp_refinement(EL);

	  if(!ER)
	    {
	      // p refinement
	      reshuffle_eta(elmt_after, nE_after, EL);
              nDOF_after += dim;
	    }
	  else
	    {
	      // h refinement
	      reshuffle_eta(elmt_after, nE_after, EL, ER);
              nDOF_after += (EL->nT - 1)*dim;
              nE_after++;
	    }
#endif
	}
    }
  else // do adaptivity according to bulk errors.
    {
      double total_err = 0;
      double bulk_err_ratio = (double)get_param("bulk_err_ratio");

      double tp;
      for(int i = 0; i < nE_before; i++)
	{
	  tp =  elmt_before[i]->eta;
	  total_err += tp*tp;
	}

#ifdef MPISRC
      tmp = 0;
      gdsum(&total_err, 1, &tmp);
#endif

      double bulk_err = 0;

      while(bulk_err < bulk_err_ratio*total_err)
	{
#ifdef MPISRC
	  tmp = DBL_MIN;
	  eta_max = elmt_after[0]->eta;
	  gdmax(&eta_max, 1, &tmp);
	  
	  if(eta_max == elmt_after[0]->eta)
	    {
	      EL = elmt_after[0];
	      
	      tp = EL->eta;
	      bulk_err += tp*tp;
	      
	      ER = hp_refinement(EL);
	      
	      if(!ER)
		{
		  // p refinement
		  reshuffle_eta(elmt_after, nE_after, EL);
		  nDOF_after += dim;
		}
	      else
		{
		  // h refinement
		  reshuffle_eta(elmt_after, nE_after, EL, ER);
		  nDOF_after += (EL->nT - 1)*dim;
		  nE_after++;
		}
	    }
	  
	  itmp = INT_MIN;
	  gimax(&nDOF_after, 1, &itmp);
	  
	  tmp = DBL_MIN;
	  gdmax(&bulk_err, 1, &tmp);
#else
	  EL = elmt_after[0];
	  
	  tp = EL->eta;
	  bulk_err += tp*tp;
	  
	  
	  ER = hp_refinement(EL);
	  
	  if(!ER)
	    {
	      // p refinement
	      reshuffle_eta(elmt_after, nE_after, EL);
	      nDOF_after += dim;
	    }
	  else
	    {
	      // h refinement
	      reshuffle_eta(elmt_after, nE_after, EL, ER);
	      nDOF_after += (EL->nT - 1)*dim;
	      nE_after++;
	    }
#endif
	}
    }


  return;
}

// deal with model approximation
void HP_adaptivity::update_Adpt_ElmtList_Model()
{ 

  if(!necessity_for_model_adpt() || !(int)get_param("is_model_adpt"))
    return;

  // first sort the elements according to theta to take care the model approximation.
  quicksort_theta_descend(elmt_after, 0, nE_after-1);

  ElementT* EL;
  ElementT* ER;

  double theta_ratio_dof = (double)get_param("theta_ratio_dof");

  int n_extra_dof = (int)ceil(theta_ratio_dof*(double)(nDOF_after-nDOF_before)) + 1;

  nDOF_model = n_extra_dof;

  printf("theta_ratio_dof is %lf, (nDOF_after-nDOF_before) is %d, # of extra dof is %d\n", theta_ratio_dof, nDOF_after-nDOF_before, n_extra_dof);

  int cnt = 0;

  int itmp;
  double tmp;
  double theta_max;

  while(cnt < n_extra_dof)
    {
#ifdef MPISRC
      theta_max = elmt_after[0]->theta;
      tmp = DBL_MIN;
      gdmax(&theta_max, 1, &tmp);

      if(theta_max == elmt_after[0]->theta)
        {
          EL = elmt_after[0];
          ER = EL->h_refinement();

          EL->theta = ER->theta = 0.0;

          reshuffle_theta(elmt_after, nE_after, EL, ER);

          ER->next = EL->next;
          EL->next = ER;

          nDOF_after += (EL->nT - 1)*dim;
          nE_after++;

          cnt +=  (EL->nT - 1)*dim;
        }

      itmp = INT_MIN;
      gimax(&nDOF_after, 1, &itmp);

      itmp = INT_MIN;
      gimax(&cnt, 1, &itmp); 
#else
      EL = elmt_after[0];
      ER = EL->h_refinement();

      EL->theta = ER->theta = 0.0;

      reshuffle_theta(elmt_after, nE_after, EL, ER);

      ER->next = EL->next;
      EL->next = ER;

      nDOF_after += (EL->nT - 1)*dim;
      nE_after++;

      cnt +=  (EL->nT - 1)*dim;
#endif
    }

  return;
}

// refine an element using h- or p-refinement.
ElementT* HP_adaptivity::hp_refinement(ElementT* EL)
{
  ElementT* ER = 0;

  int p_or_h;

  // p-refinement is allowed or not.
  int is_p_adapt = (int)get_param("is_p_adpt");

  // the limitation for the p-refinement.
  int p_max      = (int)get_param("max_p");

  // for regularity indicator alpha.
  double alpha_ratio   = (double)get_param("alpha_ratio");

  // allow multi-level h-refinement  or not.
  int i_multi_h = (int)get_param("is_multi_h_refinement");

  double alpha_threshold = (double)get_param("alpha_max");

  double alpha_min = ave < alpha_threshold? ave:alpha_threshold;

  p_or_h = 0;

  if(!is_p_adapt)
    p_or_h = 1;
  else
    {
      if(EL->alpha < alpha_ratio*alpha_min)
	{
	  if((EL->nT-1) < p_max)
	    ;
	  else
	    p_or_h = 1;
	}
      else
	p_or_h = 1;
    }

  // p refinement.
  if(p_or_h == 0)
    EL->p_refinement();

  // h refinement.
  if(p_or_h == 1)
    {
      ER = EL->h_refinement();

      ER->next = EL->next;
      EL->next = ER;
      
      if(!i_multi_h)
	EL->eta = ER->eta = 0.0;
    }

  return ER;
}

// necessity for adaptivity to deal with model approximaiton.
int HP_adaptivity::necessity_for_model_adpt()
{
  double theta_max = DBL_MIN;
  double theta_min = DBL_MAX;

  double tp_theta;

  //for(ElementT* E = elmt_after[0]; E; E = E->next)
  ElementT* E;
  for(int i = 0; i < nE_after; i++)
    {
      E = elmt_after[i];
      tp_theta  = E->theta;
      theta_max = theta_max > tp_theta? theta_max : tp_theta;
      theta_min = theta_min < tp_theta? theta_min : tp_theta;
    }

#ifdef MPISRC
  double tmp = DBL_MIN;
  gdmax(&theta_max, 1, &tmp);

  tmp = DBL_MAX;
  gdmin(&theta_min, 1, &tmp);

  ROOTONLY
    printf("theta ratio is %lf; theta_max is %lf\n", theta_max/theta_min, theta_max);
 
#endif

  double theta_threshold = (double)get_param("theta_threshold");

  int chk;

  if(theta_max/theta_min > theta_threshold)
    chk =  1;
  else
    chk = 0;

  return chk;
}

// target: p-extention and alpha. alpha measures the regularity. 
// p-extension provides an error indicator.
void HP_adaptivity::a_posteriori_error_indicators()
{
  ElementT* E;

  double* mesh = new double[nE_before+1];
  for(int i = 0; i < nE_before; i++)
    mesh[i] = elmt_before[i]->lb;
  mesh[nE_before] = elmt_before[nE_before-1]->hb;

#ifdef MPISRC
  // form the 1D smoother operator. (Only the root is used for parallel computing)
  double* Tmesh_global = 0;
  ROOTONLY
    Tmesh_global = new double[pllinfo.gne+1];

  MPI_Gatherv(mesh, nE_before, MPI_DOUBLE, Tmesh_global, pllinfo.ne, pllinfo.displs, MPI_DOUBLE, 0, 
	      MPI_COMM_WORLD);
  
  ROOTONLY
    {
      Tmesh_global[pllinfo.gne] = 1.0;
      devRecover->form_smoother_operator(pllinfo.gne, Tmesh_global);
    }
#else
  // form the 1D smoother operator.
  devRecover->form_smoother_operator(nE_before, mesh);
#endif

  // for the piecewise constant vectors.
  double** pd = dmatrix(dim, nE_before);
  double*  v  = new double[dim];

  // projection and extension from piecewise constant to linear elements
  double** p0_to_p1_prj  = dmatrix(dim, nE_before+1);

  // errors given by superconvergence and polynomial extension.
  double* err_vs_smoother  = new double[nE_before];
  double* err_vs_extension = new double[nE_before];

  // working space.
  double* work = new double[nE_before+1];

#ifdef MPISRC
  double* pd_global;
  double* work_global;
  double* prj_global;

  ROOTONLY
    {
      pd_global   = new double[pllinfo.gne];
      work_global = new double[pllinfo.gne+1];
      prj_global  = new double[pllinfo.gne+1];
    }
#endif

  double a, b, c;
  int p;
  double vp;
  double vt;
  double* vp2 = new double[2];


  // find out the minimum and maximum polynomial order.
  pmin = INT_MAX;
  for(int i = 0; i < nE_before; i++)
    {
      p = elmt_before[i]->nT - 1;
      pmin = pmin <= p? pmin : p;
    }

  // get the global minimum of p.
#ifdef MPISRC
  int itmp = INT_MAX;
  gimin(&pmin, 1, &itmp);
#endif

  pmax = INT_MIN;
  for(int i = 0; i < nE_before; i++)
    {
      p = elmt_before[i]->nT - 1;
      pmax = pmax >= p ? pmax : p;
    }

  // get the global maximum of p.
#ifdef MPISRC
  itmp = INT_MIN;
  gimax(&pmax, 1, &itmp);
#endif

  basis* btp;
  Jacobi polyL(0.0, 0.0);
  double* qpt;
  double* wgt;
  double* coef = new double[pmax+2];

  for(int pk = pmin; pk <= pmax; pk++)
    {

      // this is the p-th order derivative in each element, which is a constant vector.
      for(int i = 0; i < nE_before; i++)
	{
	  E = elmt_before[i];
	  
	  E->compt_derivative_order_p(pk, v);
      
	  for(int j = 0; j < dim; j++)
	    pd[j][i] = v[j];
	}

#ifdef MPISRC
      for(int i = 0; i < dim; i++)
	{
	  MPI_Gatherv(pd[i], nE_before, MPI_DOUBLE, pd_global, pllinfo.ne, pllinfo.displs, MPI_DOUBLE, 0, 
		      MPI_COMM_WORLD);
	  
	  // the projection and smoothing will be done at the root process. 
	  ROOTONLY
	    {
	      devRecover->smooth_p0_to_p1_projection (pd_global,   work_global);
	      devRecover->smooth_p1_to_p1_elliptic_cg(work_global, prj_global);
	    }

	  MPI_Scatterv(prj_global, pllinfo.ng, pllinfo.displs, MPI_DOUBLE, p0_to_p1_prj[i], nE_before+1, 
		       MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
#else
      // projection + smooth (through elliptic operator), based on superconvergence.
      for(int i = 0; i < dim; i++)
	{
	  devRecover->smooth_p0_to_p1_projection      (pd[i], work);
	  //devRecover->smooth_p1_to_p1_elliptic_jacobi (work,  p0_to_p1_prj[i]);
	  //devRecover->smooth_p1_to_p1_elliptic_gs     (work,  p0_to_p1_prj[i]);
	  devRecover->smooth_p1_to_p1_elliptic_cg     (work,  p0_to_p1_prj[i]);
	}
#endif
  
      // compute the error of p-th order derivative vs smoother.
      for(int i = 0; i < nE_before; i++)
	{
	  E = elmt_before[i];

	  if(pk == E->nT-1)
	    {
	      err_vs_smoother[i] = 0;
	  
	      for(int j = 0; j < dim; j++)
		{
		  //coef of linear elements.
		  a = p0_to_p1_prj[j][i];
		  b = p0_to_p1_prj[j][i+1];
	      
		  // coef of constant element
		  c = pd[j][i];
	      
		  err_vs_smoother[i] += (pow((a+b)/2-c,2.0) + pow(a-b,2.0)/12.0)*(mesh[i+1]-mesh[i]);	  
		}  
	      
	      err_vs_smoother[i] = sqrt(err_vs_smoother[i]);
	    }
	}
      
      // compute the error of p-th order derivative vs extension
      for(int i = 0; i < nE_before; i++)
	{
	  E = elmt_before[i];

	  if(pk == E->nT-1)
	    {
	      // p is the order for the polynomial extension. 
	      p = E->nT;

	      for(int j = 0; j < dim; j++)
		v[j] = (p0_to_p1_prj[j][i+1] - p0_to_p1_prj[j][i])/(mesh[i+1]-mesh[i]);
	  
	      //(p+1)-th order derivative of the mode of order p+1.
	      
	      btp = ElementT::getbasis(p+1, p+2, 't');
	      btp->getD2Cnst(vp2, mesh[i+1]-mesh[i]);
	      vp = vp2[0];

	      for(int j = 0; j < dim; j++)
		v[j] /= vp;
	      
	      // the orthogonal part of the extension to the low-order modes.
	      btp->getOrthProjection(p-1, coef);
	      dscalC(p, -1.0, coef, 1);

	      dzero(dim*(p+2), E->p_ext_GLQ[0], 1);

	      coef[p]   = coef[p-1];
	      coef[p-1] = 1.0;

	      // store the information of extension on quadrature points.
	      btp->J2Q(coef, E->p_ext_GLQ[0]);
	      for(int j = 1; j < dim; j++)
		dcopyC(p+2, E->p_ext_GLQ[0], 1, E->p_ext_GLQ[j], 1);
	 
	      for(int j = 0; j < dim; j++)
		dscalC(p+2, v[j], E->p_ext_GLQ[j], 1);

	      // record the coefficient of extention of order (p+1): this will be used by p-refinement.
	      for(int j = 0; j < dim; j++)
		dcopyC(p+1, coef, 1, E->p_extension[j], 1);
	      
	      for(int j = 0; j < dim; j++)
	        dscalC(p+1, v[j], E->p_extension[j], 1);

	      // p-th order derivative of the model of order p+1 takes the form cx, following is the 
	      // computation of the coefficient c.
	      vp = btp->getD2Linear(mesh[i+1]-mesh[i]);
	      
	      btp = ElementT::getbasis(p, p+1, 't');
	      btp->getD2Cnst(vp2, mesh[i+1]-mesh[i]);

	      if(p == 2)
		vt = coef[0]*vp2[0] + coef[p]*vp2[1];
	      else
		vt = coef[p-2]*vp2[0];
	  
	      err_vs_extension[i] = 0;
	      for(int j = 0; j < dim; j++)
		err_vs_extension[i] += (v[j]*v[j])*(vp*vp/3.0 + vt*vt)*(mesh[i+1]-mesh[i]);
	  
	      err_vs_extension[i] = sqrt(err_vs_extension[i]);
	    }
	}
    }

  // compute the scale factor alpha.
  for(int i = 0; i < nE_before; i++)
    elmt_before[i]->alpha = err_vs_smoother[i]/err_vs_extension[i];
 
  // compute the indicator theta
  compute_theta();

  // compute the indicators eta and alpha.
  compute_eta();

  delete_matrix(pd);
  delete_matrix(p0_to_p1_prj);

  delete[] v;
  delete[] work;
  delete[] err_vs_smoother;
  delete[] err_vs_extension;
  delete[] mesh;

  return;
}


// theta measures the satisfaction of arc length constraint.
void HP_adaptivity::compute_theta()
{
  ElementT* E;
  for(int i = 0; i < nE_before; i++)
    {
      E = elmt_before[i];
      E->compute_theta();
    }

  return;
}

// compute the a posteriori error indicator.
void HP_adaptivity::compute_eta()
{
  ElementT* E;
  for(int i = 0; i < nE_before; i++)
    {
      E = elmt_before[i];
      E->compute_eta();
    }

  return;
}

void HP_adaptivity::reshuffle_theta(ElementT** a, int N, ElementT* EL, ElementT* ER)
{
  for(int i = 0; i < (N-1); i++)
    a[i] = a[i+1];

  // insert EL into the list.
  int ip;
  ip = iwhere_theta_descend(a, 0, N-2, EL->theta);

  int istart;
  int iend;

  if(ip == -1)
    {
      istart = N-2;
      iend   = 0;
    }
  else
    {
      istart = N-2;
      iend   = ip+1;
    }
    
  for(int i = istart; i >= iend; i--)
    a[i+1] = a[i];

  a[iend] = EL;

  // insert ER into the list.
  ip = iwhere_theta_descend(a, 0, N-1, ER->theta);


  if(ip == -1)
    {
      istart = N-1;
      iend   = 0;
    }
  else
    {
      istart = N-1;
      iend   = ip+1;
    }


  for(int i = istart; i >= iend; i--)
    a[i+1] = a[i];

  a[iend] = ER;

  return;
}

void HP_adaptivity::reshuffle_eta(ElementT** a, int N, ElementT* EL)
{
  for(int i = 0; i < (N-1); i++)
    a[i] = a[i+1];

  a[N-1] = EL;

  return;
}

void HP_adaptivity::reshuffle_eta(ElementT** a, int N, ElementT* EL, ElementT* ER)
{
  for(int i = 0; i < (N-1); i++)
    a[i] = a[i+1];

  // insert EL into the list.
  int ip;
  ip = iwhere_eta_descend(a, 0, N-2, EL->eta);

  int istart;
  int iend;

  if(ip == -1)
    {
      istart = N-2;
      iend   = 0;
    }
  else
    {
      istart = N-2;
      iend   = ip+1;
    }

  for(int i = istart; i >= iend; i--)
    a[i+1] = a[i];

  a[iend] = EL;

  // insert ER into the list.
  ip = iwhere_eta_descend(a, 0, N-1, ER->eta);

  if(ip == -1)
    {
      istart = N-1;
      iend   = 0;
    }
  else
    {
      istart = N-1;
      iend   = ip+1;
    }
    
  for(int i = istart; i >= iend; i--)
    a[i+1] = a[i];

  a[iend] = ER;

  return;
}


void HP_adaptivity::swap_E(ElementT** left, ElementT** right)
{
  ElementT* temp;

  temp   = *left;
  *left  = *right;
  *right = temp;
}

int HP_adaptivity::partition_theta_descend(ElementT** a, int left, int right)
{
  ElementT* target = a[right];
  int i = left - 1;
  int j = right;

  while(1)
    {
      while(i<j)
	{
	  i++;
	  if(a[i]->theta <= target->theta) break;
	}

      while(j>i)
	{
	  j--;
	  if(a[j]->theta >= target->theta) break;
	}

      if(i>=j)
	break;

      swap_E(a+i, a+j);
    }

  swap_E(a+i, a+right);

  return i;
}

int HP_adaptivity::partition_eta_descend(ElementT** a, int left, int right)
{
  ElementT* target = a[right];
  int i= left - 1;
  int j= right;

  while(1)
    {
      while(i<j)
	{
          i++;
          if(a[i]->eta <= target->eta) break;
        }

      while(j>i)
	{
          j--;
	  if(a[j]->eta >= target->eta) break;
        }

      if(i>=j)
	break;

      swap_E(a+i, a+j);
    }

  swap_E(a+i, a+right);

  return i;
}


void HP_adaptivity::quicksort_theta_descend(ElementT** a, int left, int right)
{
  if(left >= right)
    return;

  int split = partition_theta_descend(a, left, right);

  quicksort_theta_descend(a, left,    split-1);
  quicksort_theta_descend(a, split+1, right  );

  return;
}

void HP_adaptivity::quicksort_eta_descend(ElementT** a, int left, int right)
{
  if(left >= right)
    return;

  int split = partition_eta_descend(a, left, right);

  quicksort_eta_descend(a, left,    split-1);
  quicksort_eta_descend(a, split+1, right  );

  return;
}



int HP_adaptivity::iwhere_theta_descend(ElementT** a, int lb, int hb, double target)
{
  int iL = lb;
  int iH = hb;
  int iM;


  if(target >= a[0]->theta)
    return -1;

  if(target <= a[hb]->theta)
    return hb;


  while(iH -iL > 1)
    {
      iM = (iH+iL)/2;
      if(target < a[iM]->theta)
        iL = iM;
      else if(target > a[iM]->theta)
        iH = iM;
      else
        return iM;
    }

  return iL;
}


int HP_adaptivity::iwhere_eta_descend(ElementT** a, int lb, int hb, double target)
{
  int iL = lb;
  int iH = hb;
  int iM;


  if(target >= a[0]->eta)
    return -1;

  if(target <= a[hb]->eta)
    return hb;

  while(iH -iL > 1)
    {
      iM = (iH+iL)/2;
      if(target < a[iM]->eta)
        iL = iM;
      else if(target > a[iM]->eta)
        iH = iM;
      else
        return iM;
    }

  return iL;
}

