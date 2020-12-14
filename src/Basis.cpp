/*
 * Basis.cpp
 *
 * Author: Xiaoliang Wan 
 */
#include <stdlib.h>
#include <math.h>
#include "FortranMapping.h"
#include "Common.h"
#include "Basis.h"

basis::basis(char bt_in, int nbm_in, int npt_in, gType gt_in)
{
  bt   = bt_in;
  nbm  = nbm_in;
  npt  = npt_in;
  gt   = gt_in;

  val  = NULL;
  D1   = NULL;
  D2   = NULL;
  Dm   = NULL;
  DmT  = NULL;
  mass = NULL;
  invm = NULL;
  wint = NULL;

  coef_orth = NULL;

  next = NULL;

  nband = -1;
}

basis::~basis(){}

void basis::setbasis()
{
  Jacobi polyL(0.0, 0.0);

  double* qpt;
  double* wgt;
  double* p1;
  double* p2;
  int i,j, k, n;
  int info;
  
  switch(bt)
    {
    case 'v':
      {
	polyL.getQW(npt, &qpt, &wgt, gt);
	
	val  = new double[nbm*npt];
	D1   = new double[nbm*npt];
	for(i = 0; i < nbm; i++)
	  getVD(npt, qpt, val+i*npt, D1+i*npt, i);
	
	Dm   = new double[npt*npt];
	DmT  = new double[npt*npt];
	polyL.getDmatrix(Dm, DmT, npt, qpt, gt);
	
	nband = 2;
	invm = new double[(nband+1)*nbm];
	dzero((nband+1)*nbm, invm, 1);
	
	for(i = 0; i < nbm; i++)
	  invm[i*(nband+1) + nband] = 2.0/(2.0*(double)i+1.0) + 
	    2.0/(2.0*(double)(i+2)+1.0);
	
	for(i = 2; i < nbm; i++)
	  invm[i*(nband+1)] = -2.0/(2.0*(double)i+1.0);
	
	dpbtrfC('U', nbm, nband, invm, nband+1, info);
	
	if(info)
	  {
	    printf("dpbtrf error when setup velocity basis!\n");
	    exit(-1);
	  }
	
	wint = new double[npt];
	dcopyC(npt, wgt, 1, wint, 1);

	break;
      }
    case 'p':
      {
	polyL.getQW(npt, &qpt, &wgt, gt);
      
	val  = new double[nbm*npt];
	D1   = new double[nbm*npt];
	for(i = 0; i < nbm; i++)
	  getVD(npt, qpt, val+i*npt, D1+i*npt, i);
	
	Dm   = new double[npt*npt];
	DmT  = new double[npt*npt];
	polyL.getDmatrix(Dm, DmT, npt, qpt, gt);
	
	nband = 0;
	invm = new double[(nband+1)*nbm];
	dzero((nband+1)*nbm, invm, 1);
	
	for(i = 0; i < nbm; i++)
	  invm[i*(nband+1) + nband] = 2.0/(2.0*(double)i+1.0);
	
	dpbtrfC('U', nbm, nband, invm, nband+1, info);
	
	if(info)
	  {
	    printf("dpbtrf error when setup velocity basis!\n");
	    exit(-1);
	  }
	
	wint = new double[npt];
	dcopyC(npt, wgt, 1, wint, 1);
	
	break;
      }
    case 't':
      /*****************************************************************
       * hp finite element basis, see Karniadakis and Sherwin's book.
       ****************************************************************/
      {
	val  = new double[nbm*npt];
	D1   = new double[nbm*npt];
	
	coef_orth = new double[nbm*nbm];
	for(i = 0; i < nbm*nbm; i++) coef_orth[i] = 0.0;
	
	mass = new double[nbm*(nbm+1)/2];
	invm = new double[nbm*(nbm+1)/2];
	
	polyL.getQW(npt, &qpt, &wgt, gt);
	
	p1 = val;
	p2 = D1;
	getVD(npt, qpt, p1, p2, 0);
	
	p1 = val + (nbm-1)*npt;
	p2 = D1  + (nbm-1)*npt;
	getVD(npt, qpt, p1, p2, nbm-1);
	
	p1 = val + npt;
	p2 = D1  + npt;
	for(i = 1; i < (nbm-1); i++)
	  {
	    getVD(npt, qpt, p1, p2, i);
	    p1 += npt;
	    p2 += npt;
	  }
	
	Dm   = new double[npt*npt];
	DmT  = new double[npt*npt];
	polyL.getDmatrix(Dm, DmT, npt, qpt, gt);
	
	for(i = 0, n = 0; i < nbm; i++)
	  for(j = i; j < nbm; j++, n++)
	    invm[n] = ddotw(npt, val+i*npt, 1, val+j*npt, 1, wgt);
	
	dcopyC(nbm*(nbm+1)/2, invm, 1, mass, 1);
	
	dpptrfC('L', nbm, invm, info);
	
	if(info)
	  {
	    printf("dpptrf error when setup time basis!\n");
	    exit(-1);
	  }
	
	wint = new double[npt];
	dcopyC(npt, wgt, 1, wint, 1);
	
	setOrthProjection();
	
	break;
      }
    default:
      ;
    }
  
  work = new double[npt];
  return;
}

// get the value and derivative of mode-th basis mode for given x.
void basis::getVD(int n, double* x, double* p, double* dp, int mode)
{
  Jacobi polyL(0.0, 0.0);
  Jacobi polyJ(1.0, 1.0);

  switch(bt)
    {
    case 'v':
      {
	double* tp  = 0;
	double* dtp = 0;
	
	if(p)
	  tp = new double[n];
	if(dp)
	  dtp = new double[n];
	
	polyL.getVD(n, x,  p,  dp, mode);
	polyL.getVD(n, x, tp, dtp, mode+2);

	if(p)
	  daxpyC(n, -1.0, tp, 1, p, 1);
	if(dp)
	  daxpyC(n, -1.0, dtp, 1, dp, 1);

	if(p)
	  delete[] tp;
	if(dp)
	  delete[] dtp;
	break;
      }
      
    case 'p':
      polyL.getVD(n, x, p, dp, mode);
      break;
      
    case 't':
      {
	double* tp = 0;

	if(!mode)
	  {
	    if(p)
	      for(int i = 0; i < n; i++)
		p[i] = (1.0 - x[i])/2.0;
	    if(dp)
	      for(int i = 0; i < n; i++)
		dp[i] = -0.5;
	  }
	else if(mode == nbm - 1)
	  {
	    if(p)
	      for(int i = 0; i < n; i++)
		p[i] = (1.0 + x[i])/2.0;
	    if(dp)
	      for(int i = 0; i < n; i++)
		dp[i] = 0.5;
	  }
	else
	  {
	    if(p && dp)
	      {
		polyJ.getVD(n, x, p, dp, mode-1);
		for(int i = 0; i < n; i++)
		  {
		    dp[i] = (1.0 - x[i]*x[i])*dp[i]/4.0 - (1.0 + x[i])*p[i]/4.0
		      + (1.0 - x[i])*p[i]/4.0;
		    p[i] *= (1.0 - x[i]*x[i])/4.0;
		  }
	      }

	    if(p && !dp)
	      {
		polyJ.getVD(n, x, p, NULL, mode-1);
		for(int i = 0; i < n; i++)
		  p[i] *= (1.0 - x[i]*x[i])/4.0;		  
	      }
	  
	    if(!p && dp)
	      {
		tp = new double[n];

		polyJ.getVD(n, x, tp, dp, mode-1);
		for(int i = 0; i < n; i++)
		  dp[i] = (1.0 - x[i]*x[i])*dp[i]/4.0 - (1.0 + x[i])*tp[i]/4.0
		    + (1.0 - x[i])*tp[i]/4.0;
		
		delete[] tp;
	      }
	  }
	break;
      }
      
    default:
      ;
    }
  
  return;
}


void basis::J2Q(double* J, double* Q)
{
  dgemvC('N', npt, nbm, 1.0, val, npt, J, 1, 0.0, Q, 1);

  return;
}


void basis::Q2J(double* Q, double* J)
{
  for(int i = 0; i < nbm; i++)
    J[i] = ddotw(npt, val+i*npt, 1, Q, 1, wint);

  int info;
  if(nband >= 0)
    dpbtrsC('U', nbm, nband, 1, invm, nband+1, J, nbm, info);
  else
    dpptrsC('L', nbm, 1, invm, J, nbm, info);

  return;
}


void basis::Grad_J2Q(double* J, double* Q)
{
  dgemvC('N', npt, nbm, 1.0, val, npt, J, 1, 0.0, work, 1);
  dgemvC('N', npt, npt, 1.0, DmT,  npt, work, 1, 0.0, Q, 1);

  return;
}

void basis::Grad_Q2Q(double* Qf, double* Qg)
{
  dgemvC('N', npt, npt, 1.0, DmT,  npt, Qf, 1, 0.0, Qg, 1);

  return;
}

// coef: projection coefficients onto previous modes.
void basis::getOrthProjection(int m, double* coef)
{
  dcopyC(m+1, coef_orth+m*nbm, 1, coef, 1);
}



void basis::setOrthProjection()
{
  if(nbm <= 2)
    return;

  double* invm_tp = new double[(nbm-1)*nbm/2];
  double* rhs = new double[nbm];

  int n = nbm*(nbm+1)/2;
  int* mask = new int[n];
  int offset;
  int ms;
  int mr;

  for(int i = 0; i < n; i++) mask[i] = 0;

  for(int k = 1; k <= nbm-2; k++)
    {
      // mark all entries related to the last mode.
      offset = nbm-1-k; mask[offset] = k;
      for(int i = 1; i < (nbm-k); i++)
	{
	  offset += (nbm-i);
	  mask[offset] = k;
	}

      for(int i = 1; i <= k; i++)
	{
	  if(!mask[offset+i])
	    mask[offset+i] = k;
	}

      ms = 0;
      mr = 0;
      for(int i = 0; i < n; i++)
	{
	  if(!mask[i])
	    invm_tp[ms++] = mass[i];
	  else
	    {
	      if(mask[i] == k)
		rhs[mr++] = mass[i];
	    }
	}

      rhs[nbm-1-k] = rhs[nbm-k];

      int info;
      dpptrfC('L', nbm-k, invm_tp, info);
      
      if(info)
	{
	  printf("dpptrf error when setup time basis!\n");
	  exit(-1);
	}
      
      dcopyC(nbm-k, rhs, 1, coef_orth+(nbm-1-k)*nbm, 1);
      dpptrsC('L', nbm-k, 1, invm_tp, coef_orth+(nbm-1-k)*nbm, nbm-k, info);
    }

  delete[] invm_tp;
  delete[] rhs;
  delete[] mask;
  return;  
}


// differentiate the last mode until a constant.
void basis::getD2Cnst(double* v, double elmt_size)
{
  int p = nbm - 1;

  if(p == 1)
    {
      v[0] = -1.0/elmt_size;
      v[1] =  1.0/elmt_size;
    }
  else
    {
      double vp = -1.0;

      for(int i = 0; i < p; i++)
	vp *= ((double)(2*p-2-i)/elmt_size);
      
      v[0] = vp;
      v[1] = 0;
    }

  return;
}


// differentiate the last mode until linear function: ax + 0.
double basis::getD2Linear(double elmt_size)
{
  int p = nbm - 2;

  double vp = -1.0;

  if(p == 1)
    vp /= elmt_size;
  else
    {
      for(int i = 0; i < p; i++)
	vp *= ((double)(2*p-i)/elmt_size);
      vp *= (double)p/2.0;
     
    }

  return vp;
}







Basis::Basis(int nX_in, int nY_in, int nptX_in, int nptY_in, char btypeY_in, gType gtypeY_in) :
  nX(nX_in), nY(nY_in), nptX(nptX_in), nptY(nptY_in), btypeY(btypeY_in), gtypeY(gtypeY_in) {}

Basis::~Basis() 
{
  delete[] bvX;
  delete[] d1X;
  delete[] d2X;

  delete[] bvY;
  delete[] d1Y;
  delete[] d2Y;
  delete[] DY;
  delete[] DtY;
  delete[] wgD;

  delete[] massY_Cholesky;

  delete[] wkF;
  delete[] wkB;
  delete[] wkD;
}

void Basis::setBasis() 
{
  int Nb;
  switch(btypeY)
    {
    case 'v':
      Nb = nY + 2;
      break;
    case 'p':
      Nb = nY;
      break;
    default:
      ;
    }

  bvX = new double[nX*nptX];
  d1X = new double[nX*nptX];
  d2X = new double[nX*nptX];
  
  bvY    = new double[Nb*nptY];
  d1Y    = new double[nY*nptY];
  d2Y    = new double[nY*nptY];
  DY     = new double[nptY*nptY];
  DtY    = new double[nptY*nptY];
  wgD    = new double[nptX*nptY];
    
  wkF = new double[nptX*nptY];
  wkB = new double[nptX*nptY];
  wkD = new double[nptX*nptY];

  if(btypeY == 'v')
    kd = 2;
  else if(btypeY == 'p')
    kd = 0;

  massY_Cholesky = new double[(kd+1)*nY];
  

  //information for the Fourier basis.  
  double x;
  for(int i = 0; i <= nX/2; i++)
    for(int j = 0; j < nptX; j++)
      {
	x = 2.0*M_PI*(double)j/(double)nptX;
	bvX[i*nptX+j] = cos((double)i*x);
      }
  for(int i = 1; i < nX/2; i++)
    for(int j = 0; j < nptX; j++)
      {
	x = 2.0*M_PI*(double)j/(double)nptX;
	bvX[(nX-i)*nptX+j] = sin((double)i*x);
      }

  for(int i = 0; i <= nX/2; i++)
    for(int j = 0; j < nptX; j++)
      {
	x = 2.0*M_PI*(double)j/(double)nptX;
	d1X[i*nptX+j] = -(double)i*sin((double)i*x);
      }
  for(int i = 1; i < nX/2; i++)
    for(int j = 0; j < nptX; j++)
      {
	x = 2.0*M_PI*(double)j/(double)nptX;
	d1X[(nX-i)*nptX+j] = (double)i*cos((double)i*x);
      }
  for(int i = 0; i <= nX/2; i++)
    for(int j = 0; j < nptX; j++)
      {
	x = 2.0*M_PI*(double)j/(double)nptX;
	d2X[i*nptX+j] = -(double)i*(double)i*cos((double)i*x);
      }
  for(int i = 1; i < nX/2; i++)
    for(int j = 0; j < nptX; j++)
      {
	x = 2.0*M_PI*(double)j/(double)nptX;
	d2X[(nX-i)*nptX+j] = -(double)i*(double)i*sin((double)i*x);
      }
  
  
  /************************************************************************* 
  Approximation space for the velocity: Due to the homogenious velocity 
  boundary conditions, the basis mode takes the combination: P_n-P_{n+2}. 
  Thus if ny is the number of basis modes in y direction, the reqiured 
  polynomial order is ny + 1.

  Approximation space for the pressure: Regular Legendre polynomials.
  *************************************************************************/
  
  Jacobi polyL(0,0);
  double* qpt;
  double* wgt;
  polyL.getQW(nptY, &qpt, &wgt, gtypeY);
  
  double* p;
  p = bvY;
  for(int i = 0; i < Nb; i++){
    polyL.getVD(nptY, qpt, p, NULL, i);
    p += nptY;
  }
  
  if(btypeY == 'v')
    {
      p = bvY;
      for(int i = 0; i < nY; i++)
	{
	  daxpyC(nptY, -1.0, p+2*nptY, 1, p, 1);
	  p += nptY;
	}
    }
  
  polyL.getDmatrix(DY, DtY, nptY, qpt, gtypeY);

  for(int i = 0; i < nY*nptY; i++) d1Y[i] = 0.0;
  for(int i = 0; i < nY*nptY; i++) d2Y[i] = 0.0;

  // compute the first-order derivative
  dgemmC('N','N', nptY, nY, nptY, 1.0, DtY, nptY, bvY, nptY, 0.0, d1Y, nptY);

  // compute the second-order derivative
  dgemmC('N','N', nptY, nY, nptY, 1.0, DtY, nptY, d1Y, nptY, 0.0, d2Y, nptY);


  // FFT from real data to Fourier coefficients in "halfcomplex" format.
  kindF = FFTW_R2HC;
  pF2B  = fftw_plan_many_r2r(1, &nptX, nptY, wkF, 0, 1, nptX, wkB, 0, 1, nptX, &kindF, FFTW_MEASURE);
  
  // inverse FFT from Fourier coefficients in "halfcomplex" format to real data.
  kindB = FFTW_HC2R;
  pB2F  = fftw_plan_many_r2r(1, &nptX, nptY, wkB, 0, 1, nptX, wkF, 0, 1, nptX, &kindB, FFTW_MEASURE);

  
  /**************************************************************
  Define the Cholesky decomposition for the mass matrix.
  For the basis P_n-P_{n+2}, the mass matrix is a 
  symmetric banded positive definite matrix with a 
  bandwidth 2.
  **************************************************************/
  int info;
      
  for(int i = 0; i < (kd+1)*nY; i++)
    massY_Cholesky[i] = 0.0;

  if(btypeY == 'v')
    {
      for(int i = 0; i < nY; i++)
	massY_Cholesky[i*(kd+1) + kd] = 2.0/(2.0*(double)i+1.0) + 
	  2.0/(2.0*(double)(i+2)+1.0);
   
      for(int i = 2; i < nY; i++)
	massY_Cholesky[i*(kd+1)] = -2.0/(2.0*(double)i+1.0);
    }

  if(btypeY == 'p')
    {
      for(int i = 0; i < nY; i++)
	massY_Cholesky[i*(kd+1) + kd] = 2.0/(2.0*(double)i+1.0);
    }

  dpbtrfC('U', nY, kd, massY_Cholesky, kd+1, info);

  // integration weights for the physical space.
  for(int i = 0; i < nptY; i++)
    for(int j = 0; j < nptX; j++)
      wgD[i*nptX+j] = 2.0*M_PI/(double)nptX*wgt[i];
  
  return;
}

// J stores the coefficients for spectral expansion, which is nY x nX.
// bvY is the Legendre basis, which is nY x nptY.
// Q stores the values on quadrature points, which is nptY x nptX. 
void Basis::J2Q(double* J, double* Q)
{  
  int NN = nptX*nptY;
  
  for(int i = 0; i < NN; i++) wkF[i] = 0.0;
  for(int i = 0; i < NN; i++) wkB[i] = 0.0;
  
  // transfer the Y direction through matrix-vector multiplication.
  dgemmC('N','T', nX, nptY, nY, 1.0, J, nX, bvY, nptY, 0.0, wkF, nX);
  
  // transfer the X direction through FFT.
  for(int i = 0; i < nptY; i++)
    cs2hc(nX, wkF+i*nX, nptX, wkB+i*nptX);
  
  // implement backward FFT transform in X direction.
  fftw_execute_r2r(pB2F, wkB, wkF);

  // copy the results to memory Q.
  dcopyC(nptX*nptY, wkF, 1, Q, 1);
}

// Q is nptY x nptX; J is nY x nX.
void Basis::Q2J(double* Q, double* J)
{
  Jacobi   lp(0,0);
  double*  npt;
  double*  wgt;
  int      info;
  
  dcopyC(nptX*nptY, Q, 1, wkF, 1);

  // implement FFT transform in X direction
  fftw_execute_r2r(pF2B, wkF, wkB);

  // transfer the hc form to cs form.
  for(int i = 0; i < nptY; i++)
    hc2cs(nptX, wkB+i*nptX, nX, wkF+i*nX);

  // get the quadrature points and integration weights.
  lp.getQW(nptY, &npt, &wgt, gtypeY);

  // projection onto each basis function in the Y direction.
  // Then solve the linear system to get the coefficients.
  for(int i = 0; i < nX; i++)
    for(int j = 0; j < nY; j++)
      *(wkB+i*nY+j) = ddotw(nptY, bvY+j*nptY, 1, wkF+i, nX, wgt);

  dpbtrsC('U', nY, kd, nX, massY_Cholesky, kd+1, wkB, nY, info);

  for(int i = 0; i < nY; i++)
    dcopyC(nX, wkB+i, nY, J+i*nX, 1);

  // scale the Fourier coefficients due to the usage of FFTW3.
  dscalC(nX*nY, 1.0/(double)nptX, J, 1);
}

void Basis::Grad_J2Q(double* J, double* dxQ, double* dyQ)
{	
  // In the Y direction, transfer it to physical space and then
  // compute the derivative using derivative matrix.
  
  // In the X direction, compute the derivative directly and then
  // implement inverse FFT.
  int NN = nptX*nptY;
  
  for(int i = 0; i < NN; i++) wkF[i] = 0.0;
  for(int i = 0; i < NN; i++) wkB[i] = 0.0;
  for(int i = 0; i < NN; i++) wkD[i] = 0.0;

  // transfer the Y direction through matrix-vector multiplication.
  dgemmC('N','T', nX, nptY, nY, 1.0, J, nX, bvY, nptY, 0.0, wkF, nX);

  if(dyQ)
    {
      // direvative in the y direction
      dgemmC('N','N', nX, nptY, nptY, 1.0, wkF, nX, DY,  nptY, 0.0, wkD, nX);
      
      // transfer the X direction through FFT to physical space.
      for(int i = 0; i < nptY; i++)
	cs2hc(nX, wkD+i*nX, nptX, wkB+i*nptX);

      // implement backward FFT in x direction
      fftw_execute_r2r(pB2F, wkB, wkD);

      dcopyC(nptX*nptY, wkD, 1, dyQ, 1);
    }

  if(dxQ)
    {
      // compute the direvative with respect to x.
      double tp;
      for(int i = 0; i < nptY; i++)
	{
	  for(int j = 1; j < nX/2; j++)
	    {
	      tp              = -wkF[i*nX+j]*(double)j;
	      wkF[i*nX+j]     = wkF[(i+1)*nX-j]*(double)j;
	      wkF[(i+1)*nX-j] = tp;
	    }
	  
	  // Derivative of mean mode is zero;
	  // Derivative of mode cos(N/2) is zero since sin(N/2) is not in the approximation space.
	  wkF[i*nX]      = 0.0;
	  wkF[i*nX+nX/2] = 0.0;
	}
      
      // transfer the X direction through FFT.
      for(int i = 0; i < nptY; i++)
	cs2hc(nX, wkF+i*nX, nptX, wkB+i*nptX);
      
      // implement backward FFT transform in X direction.
      fftw_execute_r2r(pB2F, wkB, wkF);
      
      // copy the results to memory dxQ
      dcopyC(nptX*nptY, wkF, 1, dxQ, 1);
    }
}

void Basis::Grad_Q2Q(double* Q, double* dxQ, double* dyQ)
{
  // In the Y direction, compute the derivative using the derivative matrix.
  // In the X direction, implement FFT and compute the derivative directly.
  // In the X direction, implement the inverse FFT.

  if(dxQ)
    {
      dcopyC(nptX*nptY, Q, 1, wkF, 1);

      // implement FFT transform in X direction
      fftw_execute_r2r(pF2B, wkF, wkB);
      
      // scale the Fourier coefficients due to the usage of FFTW3.
      dscalC(nptX*nptY, 1.0/(double)nptX, wkB, 1);

      // direvative in the x direction
      double tp;
      for(int i = 0; i < nptY; i++)
	{
	  for(int j = 1; j < nX/2; j++)
	    {
	      tp                 =  wkB[i*nptX+j]*(double)j;
	      wkB[i*nptX+j]      = -wkB[(i+1)*nptX-j]*(double)j;
	      wkB[(i+1)*nptX-j]  =  tp;
	    }
	  
	  // Derivative of mean mode is zero;
	  // Derivative of mode cos(N/2) is zero since sin(N/2) is not in the approximation space.
	  wkB[i*nptX]      = 0.0;
	  for(int j = nX/2; j < nptX-nX/2+1; j++)
	    wkB[i*nptX+j]  = 0.0;
	}

      fftw_execute_r2r(pB2F, wkB, wkF);

      dcopyC(nptX*nptY, wkF, 1, dxQ, 1);
    }

  if(dyQ)
    {
      // direvative in the y direction
      dgemmC('N','N', nptX, nptY, nptY, 1.0, Q, nptX, DY,  nptY, 0.0, wkF, nptX);

      dcopyC(nptX*nptY, wkF, 1, dyQ, 1);
    }
}



