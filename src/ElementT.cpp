/**********************
 * ElementT.cpp
 *
 *  Author: xlwan
 **********************/
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include "FortranMapping.h"
#include "Common.h"
#include "Basis.h"
#include "ElementT.h"

basis*  ElementT::blist = 0;

ElementT::ElementT(int nT_in, int dim_in) 
{
  nT   = nT_in;
  dim  = dim_in;

  nptT = (int)ceil((double)get_param("RNT")*(double)nT_in);
}

ElementT::~ElementT() 
{
  if(Ju) delete_matrix(Ju);
  if(Qu) delete_matrix(Qu);
  if(Du) delete_matrix(Du);

  if(fut)    delete_matrix(fut);
  if(fut_d)  delete_matrix(fut_d);
  if(fut_b)  delete_matrix(fut_b);

  if(fgt)
    {
      delete[] fgt[0][0];
      delete[] fgt[0];
      delete[] fgt;
    }

  if(loc_mapping) delete[] loc_mapping;
  if(p_extension) delete_matrix(p_extension);
  if(p_ext_GLQ)   delete_matrix(p_ext_GLQ);


  if(Imat_L) delete[] Imat_L;
  if(Imat_R) delete[] Imat_R;
}

// geometric information.
void ElementT::set_geo_info(double lb_in, double hb_in)
{
  lb = lb_in;
  hb = hb_in;
  dT = hb - lb;

  return;
}

// allocate memory for all necessary variables. 
void ElementT::mem_alloc_map()
{ 
  // Unknown coefficients.
  Ju = dmatrix(nT, dim);

  // Gradient of action functional.
  Du    = dmatrix(nT,   dim);

  // MAP in the physical space.
  Qu    = dmatrix(nptT, dim);

  // Information of the linear and nonlinear terms of dynamical system.
  fut   = dmatrix(nptT, dim);
  fut_d = dmatrix(nptT, dim);
  fut_b = dmatrix(nptT, dim);
  
  double* tp;
  double** tp2;
  tp     = new double[dim*dim*nptT];
  tp2    = new double*[dim*nptT];
  fgt    = new double**[nptT];
  fgt[0] = tp2; 
  for(int i = 0; i < nptT; i++)
    fgt[i] = fgt[0] + i*dim;

  fgt[0][0] = tp;
  for(int i = 0; i < nptT; i++)
    for(int j = 0; j < dim; j++)
      fgt[i][j] = fgt[0][0] + i*dim*dim + j*dim;

  // local mapping of Ju to a global vector used by the optimization solver.
  loc_mapping = new int[nT*dim];

  // basis with nT modes and nptT GLQ points.
  bt = ElementT::getbasis(nT, nptT, 't');

  // storage for the polynomial extension of order p+1.
  int qt_ext = nT + 2;
  p_extension = dmatrix(dim, nT+1);;
  p_ext_GLQ   = dmatrix(dim, qt_ext);

  Imat_L = new double[qt_ext*qt_ext];
  Imat_R = new double[qt_ext*qt_ext];

  Jacobi polyL(0.0, 0.0);
  double* qpt;
  double* wgt;
  polyL.getQW(qt_ext, &qpt, &wgt, Lobatto);

  // interpolation matrix to the left element.
  double* tx = new double[qt_ext];
  for(int i = 0; i < qt_ext; i++)
    tx[i] = (qpt[i]-1.0)/2.0;
  polyL.getImatrix(Imat_L, qt_ext, qpt, qt_ext, tx, Lobatto);


  // interpolation matrix to the right element.
  for(int i = 0; i < qt_ext; i++)
    tx[i] = (qpt[i]+1.0)/2.0;
  polyL.getImatrix(Imat_R, qt_ext, qpt, qt_ext, tx, Lobatto);
      
  delete[] tx;

  return;
}

void ElementT::relocate_mem_for_p_refinement()
{
  // transfer Ju and p_extension to Ju of new element.
  double** Ju_new = dmatrix(nT+1, dim);
  
  dzero((nT+1)*dim, Ju_new[0], 1);

  dcopyC((nT-1)*dim, Ju[0],    1, Ju_new[0],  1);
  dcopyC(dim,        Ju[nT-1], 1, Ju_new[nT], 1);

  // including information from p_extension may introduce instability. (not conclusive)
  /*
  for(int i = 0; i <= nT; i++)
    for(int j = 0; j < dim; j++)
      Ju_new[i][j] += p_extension[j][i];
  */

  // resize everything.
  this->~ElementT();

  nT   = nT + 1;
  nptT = (int)ceil((double)get_param("RNT")*(double)nT);
  mem_alloc_map();

  // replace Ju using Ju_new.
  delete_matrix(Ju);
  Ju = Ju_new;
  
  return;
}


// get basis; if not exists, generats it!
basis* ElementT::getbasis(int nx, int nptx, char uvp)
{
  for(basis* pb = ElementT::blist; pb; pb = pb->next)
    if(nx == pb->nbm && nptx == pb->npt && uvp == pb->bt)
      return pb;
  
  basis* pb = new basis(uvp, nx, nptx, Lobatto);
  pb->setbasis();
  
  pb->next = ElementT::blist;
  ElementT::blist = pb;
  
  return pb;
}

// prepare computation of AF.
void ElementT::pre_compt_AF()
{
  // compute the u on Quadrature points.
  for(int i = 0; i < dim; i++)
    {
      for(int j = 0; j < nptT; j++)
	{
	  Qu[j][i] = 0;
	  for(int k = 0; k < nT; k++)
	    Qu[j][i] += bt->val[k*nptT+j]*Ju[k][i];
	}
    }

  // compute the force term on Quadrature points.
  double tpx,tpy;
  double* wk = new double[dim];

  dzero(nptT*dim, fut[0],   1);
  dzero(nptT*dim, fut_d[0], 1);
  dzero(nptT*dim, fut_b[0], 1);

  double beta = 10;
  for(int i = 0; i < nptT; i++)
    {
      snapshot(i, NULL, wk);

      tpx = Qu[i][0]*Qu[i][0];
      tpy = Qu[i][1]*Qu[i][1];

      fut_d[i][0] = wk[0];
      fut_d[i][1] = wk[1];

      fut_b[i][0] = Qu[i][0]*(1-tpx-beta*tpy);
      fut_b[i][1] = -Qu[i][1]*(1+tpx);
    }

  delete[] wk;
  return;
}

// prepare computation of Gradient.
void ElementT::pre_compt_Grad()
{
  // compute the u on Quadrature points.
  for(int i = 0; i < dim; i++)
    {
      for(int j = 0; j < nptT; j++)
	{
	  Qu[j][i] = 0.0;
	  for(int k = 0; k < nT; k++)
	    Qu[j][i] += bt->val[k*nptT+j]*Ju[k][i];
	}
    }

  // compute the force term on Quadrature points.
  double tpx,tpy;
  double* wk = new double[dim];

  dzero(nptT*dim, fut[0],   1);
  dzero(nptT*dim, fut_d[0], 1);
  dzero(nptT*dim, fut_b[0], 1);

  double beta = 10;
  for(int i = 0; i < nptT; i++)
    {
      snapshot(i, NULL, wk);
      tpx = Qu[i][0]*Qu[i][0];
      tpy = Qu[i][1]*Qu[i][1];

      fut_d[i][0] = wk[0];
      fut_d[i][1] = wk[1];

      fut_b[i][0] = Qu[i][0]*(1-tpx-beta*tpy);
      fut_b[i][1] = -Qu[i][1]*(1+tpx);
    }

  double tpxy;
  for(int i = 0; i < nptT; i++)
    {
      tpx  = Qu[i][0]*Qu[i][0];
      tpy  = Qu[i][1]*Qu[i][1];
      fgt[i][0][0] = 1.0 - 3.0*tpx -beta*tpy;
    }

  for(int i = 0; i < nptT; i++)
    {
      tpxy = Qu[i][0]*Qu[i][1];
      fgt[i][0][1] = -2.0*beta*tpxy;
    }

  for(int i = 0; i < nptT; i++)
    {
      tpxy = Qu[i][0]*Qu[i][1];
      fgt[i][1][0] = -2.0*tpxy;
    }

  for(int i = 0; i < nptT; i++)
    {
      tpx  = Qu[i][0]*Qu[i][0];
      fgt[i][1][1] = -1.0-tpx;
    }

  delete[] wk;
  return;
}


// local contribution to <\dot{\phi},\dot{\phi}> - required by T(\phi).
double ElementT::local_vd()
{
  int Nt = bt->npt;
  Jacobi polyL(0.0, 0.0);
  double* qpt;
  double* wgt;

  polyL.getQW(Nt, &qpt, &wgt, Lobatto);

  double scal = dT/2.0;
  double val  = 0;

  for(int i = 0; i < Nt; i++)
    val += (fut_d[i][0]*fut_d[i][0]+fut_d[i][1]*fut_d[i][1])*wgt[i]*scal;

  return val;
}

// local contribution to <b(\phi),b(\bphi)> - required by T(\phi).
double ElementT::local_vb()
{
  int Nt = bt->npt;
  Jacobi polyL(0.0, 0.0);
  double* qpt;
  double* wgt;
  
  polyL.getQW(Nt, &qpt, &wgt, Lobatto);
  
  double scal = dT/2.0;
  double val  = 0;

  for(int i = 0; i < Nt; i++)
    val += (fut_b[i][0]*fut_b[i][0]+fut_b[i][1]*fut_b[i][1])*wgt[i]*scal;

  return val;
}

// local contribution to Action Functional.
double ElementT::local_val_AF()
{
  int Nt = bt->npt;
  Jacobi polyL(0.0, 0.0);
  double* qpt;
  double* wgt;

  polyL.getQW(Nt, &qpt, &wgt, Lobatto);

  double scal = dT/2.0;
  double val  = 0;

  double tp0, tp1;

  for(int i = 0; i < Nt; i++)
    {
      tp0 = fut_d[i][0]/T - fut_b[i][0];
      tp1 = fut_d[i][1]/T - fut_b[i][1];
      val += (tp0*tp0 + tp1*tp1)*wgt[i]*scal;

    }
  return val/2.0*T;
}

// local contribution to gradient of action funcitonal.
void ElementT::local_grad_AF()
{
  for(int i = 0; i < nT; i++)
    for(int j = 0; j < dim; j++)
      Du[i][j] = delta_grad(i,j);
}

// indicator of the distant from arc length constraint. 
double ElementT::local_arc_length_constraint()
{
  Jacobi polyL(0.0, 0.0);
  double* qpt;
  double* wgt;

  polyL.getQW(nptT, &qpt, &wgt, Lobatto);

  double scal = dT/2.0;
  double val  = 0;

  double tp0, tp1;

  for(int i = 0; i < nptT; i++)
    {
      tp0 = sqrt((fut_d[i][0]*fut_d[i][0] + fut_d[i][1]*fut_d[i][1])/T/T);
      tp1 = sqrt(fut_b[i][0]*fut_b[i][0] + fut_b[i][1]*fut_b[i][1]);
      val += (tp0-tp1)*(tp0-tp1)*wgt[i]*scal;
    }

  return sqrt(val*T);
}

// gradient of a certain mode.
double ElementT::delta_grad(int iT, int jU)
{
  int Nt = bt->npt;

  double* qpt;
  double* wgt;
  Jacobi polyL(0.0, 0.0);  

  polyL.getQW(Nt, &qpt, &wgt, Lobatto);

  double val = 0;
  
  if(jU == 0)
    for(int i = 0; i < nptT; i++)
      {
	val += T*(fut_d[i][0]/T - fut_b[i][0])*(bt->D1[iT*nptT+i]*2.0/dT/T 
						- fgt[i][0][0]*bt->val[iT*nptT+i])*wgt[i]*dT/2.0;
	val += T*(fut_d[i][1]/T - fut_b[i][1])*(-fgt[i][1][0]*bt->val[iT*nptT+i])*wgt[i]*dT/2.0;
      }

  if(jU == 1)
    for(int i = 0; i < nptT; i++)
      {
	val += T*(fut_d[i][0]/T - fut_b[i][0])*(-fgt[i][0][1]*bt->val[iT*nptT+i])*wgt[i]*dT/2.0;
	val += T*(fut_d[i][1]/T - fut_b[i][1])*(bt->D1[iT*nptT+i]*2.0/dT/T 
						- fgt[i][1][1]*bt->val[iT*nptT+i])*wgt[i]*dT/2.0;
      }
  
  return val;
}


// compute arc length.
double ElementT::get_arc_length(double t)
{
  if(lb - t > 10e-6 || t - hb > 10e-6)
    {
      printf("get_arc_length: t is out of range!\n");
      exit(-1);
    }

  int Nt = bt->npt;

  double* qpt;
  double* wgt;
  Jacobi polyL(0.0, 0.0);
  polyL.getQW(Nt, &qpt, &wgt, Lobatto);
 
  double* vt = new double[dim];
  double tp;

  double val = 0;
  for(int i = 0; i < Nt; i++)
    {
      tp = (t-lb)/2.0*qpt[i] + (t+lb)/2.0;
      
      snapshot(tp, NULL, vt);

      val += sqrt(vt[0]*vt[0] + vt[1]*vt[1])*wgt[i]*(t-lb)/2.0;
    }

  delete[] vt;

  return val;
} 


// split the current element [lb hb] to two equidistant elements.
// The current element will be [lb (lb+hb)/2]; the newly generated 
// element will be [(hb+lb)/2 hb].  
ElementT* ElementT::h_refinement()
{
  /*************** First deal with the new element (Right)*****************/
  ElementT* new_elmt = new ElementT(nT, dim);

  // allocate memory
  new_elmt->mem_alloc_map();

  // set geometric information
  new_elmt->set_geo_info((hb+lb)/2, hb);

  // get the Gauss-Lobatto quadrature points.
  double* qpt;
  double* wgt;
  Jacobi polyL(0.0, 0.0);
  polyL.getQW(nptT, &qpt, &wgt, Lobatto);
 
  // projection the current solution onto the new element.
  double tp;
  double a = (hb+lb)/2.0;
  double b = hb;

  for(int i = 0; i < nptT; i++)
    {
      tp = (b-a)/2.0*qpt[i] + (b+a)/2.0;
      snapshot(tp, new_elmt->Qu[i], NULL);
    }
  new_elmt->Q2J();

  // interpolate the extention of (p+1)-th order onto the new element.
  int qt_ext = nT + 2;
  for(int i = 0; i < dim; i++)
    dgemvC('T', qt_ext, qt_ext, 1.0, Imat_R, qt_ext, p_ext_GLQ[i], 1, 0.0, new_elmt->p_ext_GLQ[i], 1); 


  // information inherited from parent element.
  new_elmt->alpha = alpha;
  new_elmt->T     = T;
  
  new_elmt->pre_compt_AF();
  new_elmt->compute_eta();
  new_elmt->compute_theta();

  /***************** Then deal with the old element (Left) ***********************/
  a = lb;
  b = (hb+lb)/2.0;

  for(int i = 0; i < nptT; i++)
    {
      tp = (b-a)/2.0*qpt[i] + (b+a)/2.0;
      snapshot(tp, Qu[i], NULL);
    }
 
  Q2J();
  
  double* wk = new double[qt_ext];
  for(int i = 0; i < dim; i++)
    {
      dgemvC('T', qt_ext, qt_ext, 1.0, Imat_L, qt_ext, p_ext_GLQ[i], 1, 0.0, wk, 1);
      dcopyC(qt_ext, wk, 1, p_ext_GLQ[i], 1);
    }
  delete[] wk;

  //update geometric information
  hb = (hb+lb)/2.0;
  dT = hb - lb;
  
  pre_compt_AF();
  compute_eta();
  compute_theta(); 

  return new_elmt;
}

void ElementT::p_refinement()
{
  // relocate the memeory
  relocate_mem_for_p_refinement();

  // set the error indicator to be 0, i.e., p-refinement is vaid for one level.
  eta = 0.0;

  return;
}

// compute the indicator for the linear time scaling
void ElementT::compute_theta()
{
  theta = local_arc_length_constraint();

  return;
}

// compute the error indicator eta assuming that alpha and p_extension are known.
void ElementT::compute_eta()
{
  int qt_ext  = nT + 2; 

  basis* btp = ElementT::getbasis(nT+1, qt_ext, 't');

  double* wk = new double[qt_ext];

  // get the Gauss-Lobatto quadrature points.
  double* qpt;
  double* wgt;
  Jacobi polyL(0.0, 0.0);
  polyL.getQW(qt_ext, &qpt, &wgt, Lobatto);
  

  // alpha * |phi_{p+1}-\phi_{p}|_1
  eta = 0;
  for(int i = 0; i < dim; i++)
    {
      btp->Grad_Q2Q(p_ext_GLQ[i], wk);
      eta += ddotw(qt_ext, wk, 1,wk, 1, wgt)*(2.0/dT);
    }
  
  eta  = sqrt(eta);
  eta *= alpha;

  return;
}


// compute the derivative of the solution of order nT-1. 
void ElementT::compt_derivative_order_p(int p, double* vc)
{
  if( p > nT-1)
    {
      dzero(dim, vc, 1);
      return;
    }

  // only linear modes.
  if(p == 1)
    {
      for(int i = 0; i < dim; i++)
	vc[i] = (Ju[nT-1][i]-Ju[0][i])/dT;
    }
  else
    {
      // high-order modes with highest order nT-1.
      for(int i = 0; i < dim; i++)
	{
	  vc[i] = Ju[p-1][i];
	  
	  for(int j = 0; j < p; j++)
	    vc[i] *= ((double)(2*p-2-j)/dT);
	
	  vc[i] *= -1.0;
	}
    }

  return;
}




// modal space to physical space.
void ElementT::J2Q()
{
  // compute u on Quadrature points.
  double* Jtp = new double[nT];
  double* Qtp = new double[nptT];

  for(int i = 0; i < dim; i++)
    {
      dcopyC(nT, Ju[0]+i, dim, Jtp, 1);
      bt->J2Q(Jtp, Qtp);
      dcopyC(nptT, Qtp, 1, Qu[0]+i, dim);
    }


  delete[] Jtp;
  delete[] Qtp;
  return;
}


// physical space to modal space. 
void ElementT::Q2J()
{
  double* Jtp = new double[nT];
  double* Qtp = new double[nptT];

  for(int i = 0; i < dim; i++)
    {
      dcopyC(nptT, Qu[0]+i, dim, Qtp, 1);
      bt->Q2J(Qtp, Jtp);
      dcopyC(nT, Jtp, 1, Ju[0]+i, dim);
    }

  delete[] Jtp;
  delete[] Qtp;
  return;
}


//snapshot of the solution for a certain quadrature point.
void ElementT::snapshot(int in, double* p, double* pt)
{
  int Nt = bt->npt;

  if(p)
    {
      dzero(dim, p, 1);
      for(int i = 0; i < nT; i++)
	daxpyC(dim, bt->val[i*Nt+in], Ju[i], 1, p, 1);
    }

  if(pt)
    {
      dzero(dim, pt, 1);
      for(int i = 0; i < nT; i++)
	daxpyC(dim, bt->D1[i*Nt+in]*2.0/dT, Ju[i], 1, pt, 1);
    }
  return;
}

// snapshot of the solution at a arbitrarily given time. 
void ElementT::snapshot(double t, double* v, double* vt)
{
  if(v)
    dzero(dim, v, 1);

  if(vt)
    dzero(dim, vt, 1);

  double tp  = 2.0*t/(hb-lb) - (hb+lb)/(hb-lb);
  double tv  = 0;
  double tvt = 0;

  for(int i = 0; i < nT; i++)
    {
      bt->getVD(1, &tp, &tv, &tvt, i);
      if(v)
	daxpyC(dim, tv,  Ju[i], 1, v,  1);
      if(vt)
	daxpyC(dim, tvt*2.0/dT, Ju[i], 1, vt, 1);
    }

  return;
}

// transfer a given function to Ju.
void ElementT::func2Ju(double (**f)(double))
{
  Jacobi polyL(0.0, 0.0);
  double* qpt;
  double* wgt;

  polyL.getQW(nptT, &qpt, &wgt, Lobatto);

  double tp;
  for(int i = 0; i < nptT; i++)
    {
      tp = (hb - lb)/2.0*qpt[i] + (hb + lb)/2.0;
      for(int j = 0; j < dim; j++)
        {
	  Qu[i][j] = f[j](tp);
        }
    }
  
  double* v = new double[nT];

  int info;
  for(int i = 0; i < dim; i++)
    {
      dzero(nT, v, 1);
      for(int j = 0; j < nT; j++)
	for(int k = 0; k < nptT; k++)
	    v[j] += Qu[k][i]*bt->val[j*nptT+k]*wgt[k];
      dpptrsC('L', nT, 1, bt->invm, v, nT, info);
      dcopyC(nT, v, 1, Ju[0]+i, dim);
    }

  delete[] v;
  return;
}


