/*
 * Smoother1D.cpp
 *
 * MovingMesh.cpp provides a second-order finite element solver for 
 * the elliptic equation to adjust the mesh.
 *
 *  Author: Xiaoliang Wan
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "FortranMapping.h"
#include "Common.h"
#include "Smoother1D.h"

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
  if(xb) delete[] xb;
  if(rb) delete[] rb;
  if(Q)  delete_matrix(Q);
  if(A)  delete_matrix(A);
}

void Smoother1D::form_smoother_operator(int nTE_in, double* mesh)
{
  nTE = nTE_in;


  int n = nTE+1;
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
      dT = new double[nTE];
    }
  else
    dT = new double[nTE];

  for(int i = 0; i < nTE; i++)
    dT[i] = mesh[i+1] - mesh[i];

  dzero(3*n, Q[0], 1);
  dzero(3*n, A[0], 1);

  // mass matrix - tridiagonal.
  for(int i = 1; i < n; i++)
    Q[0][i] = dT[i-1]/6.0;

  Q[1][0] = dT[0]/3.0;
  for(int i = 1; i < n-1; i++)
    Q[1][i] = (dT[i-1]+dT[i])/3.0;
  Q[1][n-1] = dT[nTE-1]/3.0;

  for(int i = 0; i < (n-1); i++)
    Q[2][i] = dT[i]/6.0;
  
  // stiffness matrix - tridiagonal.
  for(int i = 1; i < n; i++)
    A[0][i] = -1.0/dT[i-1] + Q[0][i];

  A[1][0] = 1.0/dT[0] + Q[1][0];
  for(int i = 1; i < n-1; i++)
    A[1][i] = (1.0/dT[i-1] + 1.0/dT[i]) + Q[1][i];
  A[1][n-1] = 1.0/dT[nTE-1] + Q[1][n-1];

  for(int i = 0; i < n-1; i++)
    A[2][i] = -1.0/dT[i] + Q[2][i];

  return;
}

//smooth piecewise constant - through projection
void Smoother1D::smooth_p0_to_p1_projection(double* p0, double* p1)
{
  int n = nTE+1;

  rb[0] = p0[0]*dT[0]/2.0;
  for(int i = 1; i < n-1; i++)
    rb[i] = (p0[i-1]*dT[i-1] + p0[i]*dT[i])/2.0;
  rb[n-1] = p0[nTE-1]*dT[nTE-1]/2.0;

  tridag_solver(Q[0], Q[1], Q[2], rb, p1, n);

  return;
}

//smooth piecewise linear element using conjugate gradient
void Smoother1D::smooth_p1_to_p1_elliptic_cg(double* p1_before, double* p1_after)
{
  int n = nTE+1;

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
  int n = nTE+1;

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
  int n = nTE+1;

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
  int n = nTE+1;

  x_out[0] = A[1][0]*x_in[0] + A[0][1]*x_in[1];
  for(int i = 1; i < n-1; i++)
    x_out[i] = A[2][i-1]*x_in[i-1]+A[1][i]*x_in[i]+A[0][i+1]*x_in[i+1];
  x_out[n-1]= A[2][n-2]*x_in[n-2] + A[1][n-1]*x_in[n-1];

  return;
}

// matrix vector multiplication: Rx - R is A without the diagonal line.
void Smoother1D::Rxv(double* x_in, double* x_out)
{
  int n = nTE+1;

  x_out[0] = A[0][1]*x_in[1];
  for(int i = 1; i < n-1; i++)
    x_out[i] = A[2][i-1]*x_in[i-1]+A[0][i+1]*x_in[i+1];
  x_out[n-1]= A[2][n-2]*x_in[n-2];

  return;
}
