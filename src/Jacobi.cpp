/*
 * Legendre.cpp
 *
 *  Created on: May 31, 2010
 *      Author: xlwan
 */
#include <stdlib.h>
#include <stdio.h>
#include "FortranMapping.h"
#include "Jacobi.h"
#include "polylib.h"

Qstore* Jacobi::qw = 0;

Jacobi::Jacobi(double alpha_in, double beta_in) 
{
  alpha = alpha_in;
  beta  = beta_in;
}

Jacobi::~Jacobi() {}

void Jacobi::getQW(int npt, double** qpt, double** wgt, gType type) 
{
  
  for(Qstore* bp = Jacobi::qw; bp; bp = bp->next)
    {
      if(bp->alpha == alpha && bp->beta == beta && bp->n == npt && bp->type == type)
	{
	  *qpt = bp->qpt;
	  *wgt = bp->wgt;

	  return;
	}
    }
  
  Qstore* bp = (Qstore*)malloc(sizeof(Qstore));
  bp->alpha = alpha;
  bp->beta  = beta;
  bp->n     = npt;
  bp->qpt   = new double[npt];
  bp->wgt   = new double[npt];
  bp->type  = type;
  bp->next  = Jacobi::qw;
  
  Jacobi::qw = bp;
  
  switch (type) 
    {
    case Gauss: // Gauss.
      zwgj(bp->qpt, bp->wgt, npt, alpha, beta);
      break;
    case Lobatto: // Gauss-Lobatto.
      zwglj(bp->qpt, bp->wgt, npt, alpha, beta);
      break;
    case RadauL: // Gauss Radau with left end.
      zwgrjm(bp->qpt, bp->wgt, npt, alpha, beta);
      break;
    case RadauR: // Gauss Radau with right end.
      zwgrjp(bp->qpt, bp->wgt, npt, alpha, beta);
      break;
    }
  
  *qpt = bp->qpt;
  *wgt = bp->wgt;
 
  return;
}

void Jacobi::getVD(int nx, double* x, double* val, double* div, int n) 
{
  if (val)
    jacobfd(nx, x, val, NULL, n, alpha, beta);

  if(div)
    jacobd(nx, x, div, n, alpha, beta);
  
  return;
}

void Jacobi::getDmatrix(double* D, double* Dt, int npt, double* qpt, gType type) 
{
  switch (type) 
    {
    case Gauss:
      Dgj(D, Dt, qpt, npt, alpha, beta);
      break;
    case Lobatto:
      Dglj(D, Dt, qpt, npt, alpha, beta);
      break;
    case RadauL:
      Dgrjm(D, Dt, qpt, npt, alpha, beta);
      break;
    case RadauR:
      Dgrjp(D, Dt, qpt, npt, alpha, beta);
      break;
    }
  
  return;
}

void Jacobi::getImatrix(double *Q2X, int npt, double *qpt, int nx, double *x, gType type) 
{
  switch (type) 
    {
    case Gauss:
      Imgj(Q2X, qpt, x, npt, nx, alpha, beta);
      break;
    case Lobatto:
      Imglj(Q2X, qpt, x, npt, nx, alpha, beta);
      break;
    case RadauL:
      Imgrjm(Q2X, qpt, x, npt, nx, alpha, beta);
      break;
    case RadauR:
      Imgrjp(Q2X, qpt, x, npt, nx, alpha, beta);
      break;
    }
  
  return;
}

