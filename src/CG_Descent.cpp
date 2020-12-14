/*************************************************************************** 
   =========================================================================
   ============================ CG_DESCENT =================================
   =========================================================================
       ________________________________________________________________
      |      A conjugate gradient method with guaranteed descent       |
      |             C-code Version 1.1  (October 6, 2005)              |
      |                    Version 1.2  (November 14, 2005)            |
      |                    Version 2.0  (September 23, 2007)           |
      |                    Version 3.0  (May 18, 2008)                 |
      |           William W. Hager    and   Hongchao Zhang             |
      |          hager@math.ufl.edu       hzhang@math.ufl.edu          |
      |                   Department of Mathematics                    |
      |                     University of Florida                      |
      |                 Gainesville, Florida 32611 USA                 |
      |                      352-392-0281 x 244                        |
      |                                                                |
      |                 Copyright by William W. Hager                  |
      |                                                                |
      |          http://www.math.ufl.edu/~hager/papers/CG              |
      |________________________________________________________________|
       ________________________________________________________________
      |This program is free software; you can redistribute it and/or   |
      |modify it under the terms of the GNU General Public License as  |
      |published by the Free Software Foundation; either version 2 of  |
      |the License, or (at your option) any later version.             |
      |This program is distributed in the hope that it will be useful, |
      |but WITHOUT ANY WARRANTY; without even the implied warranty of  |
      |MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the   |
      |GNU General Public License for more details.                    |
      |                                                                |
      |You should have received a copy of the GNU General Public       |
      |License along with this program; if not, write to the Free      |
      |Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, |
      |MA  02110-1301  USA                                             |
      |________________________________________________________________|

  >>  Modified by Xiaoliang Wan.

***************************************************************************/

#include "CG_Descent.h"
#include "Common.h"


CG_Descent::CG_Descent(){}

CG_Descent::CG_Descent(int n)
{
  _work = (double*)malloc(9*n*sizeof(double)) ;

  _Com = (cg_com*)malloc(sizeof(cg_com)) ;

  _Com->d       = _work + 2*n ;
  _Com->g       = _work + 3*n ;
  _Com->xtemp   = _work + 4*n ;
  _Com->gtemp   = _work + 5*n ;
  _Com->d_invP  = _work + 6*n ;
  _Com->g_P     = _work + 7*n ;
  _Com->gtemp_P = _work + 8*n ;

  _Com->n     = n ;
  _Com->nf    = (int)0 ;
  _Com->ng    = (int)0 ;
}


void CG_Descent::cg_locate_mem(int n)
{
  if(_work) free(_work);
  if(_Com)  free(_Com);

  _work = (double*)malloc(9*n*sizeof(double)) ;
  _Com  = (cg_com*)malloc(sizeof(cg_com)) ;

  _Com->d       = _work + 2*n ;
  _Com->g       = _work + 3*n ;
  _Com->xtemp   = _work + 4*n ;
  _Com->gtemp   = _work + 5*n ;
  _Com->d_invP  = _work + 6*n ;
  _Com->g_P     = _work + 7*n ;
  _Com->gtemp_P = _work + 8*n ;

  _Com->n   = n ;
  _Com->nf  = (int)0 ;
  _Com->ng  = (int)0 ;

  return;
}


int CG_Descent::cg_set_parameter
(
 cg_parameter* UParm,
 void*      instance,
 double   (*value) (void*, double*, int),
 void      (*grad) (void*, double*, double*, int),
 double (*valgrad) (void*, double*, double*, int),
 void   (*preconU) (void *, double *, double *, int, int)
 )
{
  if ( !UParm )
    {
      _Com->Parm = (cg_parameter*)malloc(sizeof(cg_parameter)) ;
      cg_default(_Com->Parm) ;
    }
  else
    _Com->Parm = UParm ;

  _Com->AWolfe     = _Com->Parm->AWolfe ;
  _Com->instance   = instance ;
  _Com->cg_value   = value ;
  _Com->cg_grad    = grad ;
  _Com->cg_valgrad = valgrad ;
  _Com->cg_preconU = preconU;

  return 1;
}

int CG_Descent::cg_initial_work(double* x, double grad_tol, int* maxit, int* status_code)
{  
  int i ;
  int n = _Com->n ;

  _Com->x = x ;

  _StopRule = _Com->Parm->StopRule ;
 
  // To synchronize each cpu, we need to set the same restart number; 
  // otherwise, parallel code may stall when each cpu has different 
  // number of unknowns.
#ifdef MPISRC
  _nrestart = (int)((double)(pllinfo.nd*pllinfo.gne)*_Com->Parm->restart_fac) ;
#else
  _nrestart = (int)(((double)n)*_Com->Parm->restart_fac) ;
#endif

  if(_Com->Parm->maxit_fac == INF)
    *maxit = _maxit = INT_INF ;
  else
    *maxit = _maxit = (int)(((double)n)*_Com->Parm->maxit_fac) ;

  _f  = ZERO ;
  _Ck = ZERO ;
  _Qk = ZERO ;


  double*   g = _Com->g ;
  double* g_P = _Com->g_P ;
  double*   d = _Com->d ;

  double xnorm = cg_linf(x, n) ;

  _f        = cg_fg (g, x, _Com) ;
  _Com->f0  = 2.0*_f ;

  // preconditioned gradient.
  cg_precon ( g, g_P, _Com ) ;  

  _gnorm    = cg_linf ( g,   n ) ;
  _gnorm_P  = cg_linf ( g_P, n ) ;

  _gnorm2   = cg_dot ( g,   g,   n ) ;
  _gnorm2_P = cg_dot ( g,   g_P, n ) ;

  // inital direction.
  for ( i = 0; i < n; i++ )
    d[i] = -g_P[i] ;

  /* check that starting function value is nan */
  if ( _f != _f )
    {
      *status_code = -1 ;
      return 0 ;
    }

  _grad_tol = grad_tol ;

  if ( _Com->Parm->StopRule ) 
    _tol = MAX (_gnorm*_Com->Parm->StopFac, grad_tol) ;
  else                  
    _tol = grad_tol ;

  if ( _Com->Parm->PrintLevel >= 1 )
#ifdef MPISRC
    ROOTONLY   
#endif
      printf ( "iter: %5i f = %14.6e gnorm = %14.6e AWolfe = %2i\n",
	       (int) 0, _f, _gnorm, _Com->AWolfe ) ; 

  if ( cg_tol (_f, _gnorm, _StopRule, _tol) )
    {
      *status_code = 0 ;
      return 0 ;
    }

  _dphi0  = -_gnorm2_P ;
  _delta2 = 2.0*_Com->Parm->delta - ONE ;
  _eta_sq = _Com->Parm->eta*_Com->Parm->eta ;
  _alpha  = _Com->Parm->step ;

  if ( _alpha == 0. )
    {
      _alpha = _Com->Parm->psi0*xnorm/_gnorm_P ;

      if ( xnorm == ZERO )
        {
	  if ( _f != ZERO ) 
	    {
	      _alpha = _Com->Parm->psi0*fabs (_f)/_gnorm2_P ; 
	    }
	  else             
	    _alpha = ONE ;
        }
    }

  return 1;
}

int CG_Descent::cg_forward_one(double* x, int* status_code)
{
  int i,status ;
  double t ;
  cg_parameter *Parm = _Com->Parm ;
  double talpha,ftemp,denom ;
  double* xtemp = _Com->xtemp ;
  double* d     = _Com->d ;
  int n         = _Com->n ;

  _Com->QuadOK  = FALSE ;
  _alpha = Parm->psi2*_alpha ;

  if ( Parm->QuadStep )
    {
      if( _f != ZERO ) 
	t = fabs ((_f-_Com->f0)/_f ) ;
      else             
	t = ONE ;
     
      if( t > Parm->QuadCutOff )       /* take provisional step talpha */
	{
	  talpha = Parm->psi1*_alpha ;
	  cg_step (xtemp, x, d, talpha, n) ;
	  ftemp = cg_f (xtemp, _Com) ;  /* provisional function value */

	  /* check if function value is nan */
	  if ( ftemp != ftemp ) /* reduce stepsize */
	    {
	      for ( i = 0; i < Parm->nexpand; i++ )
		{
		  talpha /= Parm->rho ;
		  cg_step ( xtemp, x, d, talpha, n) ;
		  ftemp = cg_f (xtemp, _Com) ;
		 
		  if ( ftemp == ftemp ) 
		    break ;
		}
	      if ( i == Parm->nexpand )
		{
		  *status_code = -2 ;
		  return 0 ;
		}
	    }

	  if ( ftemp < _f )              /* check if quadstep > 0 */
	    {
	      denom = 2.*(((ftemp-_f)/talpha)-_dphi0) ;
	      if ( denom > ZERO )        /* try a quadratic fit step */
		{
		  _Com->QuadOK = TRUE ;
		  _alpha = -_dphi0*talpha/denom ;
		}
	    }
	}
    }
  
  _Com->f0 = _f ;                          /* f0 saved as prior value */

  if ( Parm->PrintLevel >= 1 )
#ifdef MPISRC
    ROOTONLY
#endif    
      printf ("QuadOK: %2i initial a: %14.6e f0: %14.6e dphi: %14.6e\n",
	      _Com->QuadOK, _alpha, _Com->f0, _dphi0) ;
    
  /* parameters in Wolfe and approximate Wolfe conditions, and in update */

  _Qk = Parm->Qdecay*_Qk + ONE ;
  _Ck = _Ck + (fabs (_f) - _Ck)/_Qk ;        /* average cost magnitude */

  if ( Parm->PertRule ) 
    _Com->fpert = _f + Parm->eps*_Ck ;
  else                
    _Com->fpert = _f + Parm->eps ;

  _Com->wolfe_hi = Parm->delta*_dphi0 ;
  _Com->wolfe_lo = Parm->sigma*_dphi0 ;
  _Com->awolfe_hi = _delta2*_dphi0 ;
  _Com->alpha    = _alpha ;        /* either prior step or quadratic fit step */
  _Com->f = _f ;

  if ( _Com->AWolfe ) 
    status = cg_line (_dphi0, _Com) ; /* approx. Wolfe */
  else              
    status = cg_lineW (_dphi0, _Com) ;/* ordinary Wolfe */

  if ( (status > 0) && !_Com->AWolfe )/*try approximate Wolfe line search*/
    {
      if ( Parm->PrintLevel >= 1 )
#ifdef MPISRC
	ROOTONLY
#endif
	  printf ("\nWOLFE LINE SEARCH FAILS\n") ;

      _Com->AWolfe = TRUE ;

      status = cg_line (_dphi0, _Com) ;
    }

  _alpha = _Com->alpha ;
  _f     = _Com->f ;
  _dphi  = _Com->df ;

  if ( status ) 
    {
      *status_code = status;
      return 0;
    }
  
  /*Test for convergence to within machine epsilon    
    [set feps to zero to remove this test] */
  if ( -_alpha*_dphi0 <= Parm->feps*fabs (_f) )
    {
      *status_code = 1 ;
      return 0;
    }

  return 1;
}


int CG_Descent::cg_update_direction(double* x, int iter, int* status_code)
{
  int i;
  double t,gnorm,ykyk,ykgk,dkyk;
  double beta = 0;

  double dnorm2_P;
  
  int nrestart = _nrestart ;

  int n = _Com->n ;
  double* xtemp = _Com->xtemp ;
  double* d     = _Com->d ;
  double* gtemp = _Com->gtemp ;
  double* g     = _Com->g ;
  cg_parameter* Parm = _Com->Parm ;

  double* d_invP = _Com->d_invP ; 
  double* gtemp_P= _Com->gtemp_P ;
  double* g_P    = _Com->g_P ;

  double* wk_x = _work ;
  double* wk_y = _work + n ;


  if ( iter % nrestart != 0 )
    {
      cg_precon_inv ( d, d_invP, _Com ) ;
      cg_precon ( gtemp, gtemp_P, _Com ) ;
      cg_precon ( g,     g_P,     _Com ) ;

      cg_copy (x, xtemp, n) ;

      dnorm2_P = cg_dot(d, d_invP, n);
      _gnorm   = cg_linf(gtemp, n);

      cg_dvsub ( n, gtemp,   g,   wk_x) ;
      cg_dvsub ( n, gtemp_P, g_P, wk_y) ;

      ykgk = cg_dot ( wk_x, gtemp_P, n ) ;
      ykyk = cg_dot ( wk_x, wk_y,    n ) ;
	
      // copy gtemp to g.
      cg_copy ( g, gtemp, n ) ;

      if ( cg_tol (_f, _gnorm, _StopRule, _tol) )
	{
	  *status_code = 0 ;
	  return 0 ;
	}

      dkyk = _dphi - _dphi0 ;
      //beta = (ykgk - 2.*_dphi*ykyk/dkyk)/dkyk ;
      beta = (ykgk - _dphi*ykyk/dkyk)/dkyk ;

      /* 
	 faster: initialize dnorm2 = gnorm2 at start, then      
	 dnorm2 = gnorm2 + beta**2*dnorm2 - 2.*beta*dphi                     
	 gnorm2 = ||g_{k+1}||^2                                                       
	 dnorm2 = ||d_{k+1}||^2                                                       
	 dpi = g_{k+1}' d_k 
      */      
      t = -ONE/sqrt (dnorm2_P*MIN (_eta_sq, _gnorm2_P)) ;
      
      beta = MAX (beta, t) ;

      /*    update search direction d = -g + beta*dold */
      cg_precon ( g, g_P, _Com ) ;
      
      _gnorm2_P = cg_dot ( g, g_P, n ) ;

      for(i = 0; i < n; i++)
	d[i] = -g_P[i] + beta*d [i] ;
      
      _dphi0 = -_gnorm2_P + beta*_dphi ;

      if ( Parm->debug ) /* Check the dphi0 = d'g */
	{
	  t = ZERO ;
	  for (i = 0; i < n; i++)  
	    t = t + d [i]*g [i] ;

	  if ( fabs(t-_dphi0) > Parm->debugtol*fabs(_dphi0) )
	    {
	      printf("Warning, dphi0 != d'g!\n");
	      printf("dphi0:%14.6e, d'g:%14.6e\n",_dphi0, t) ;
	    }
	}
    }
  else
    {
      /* search direction d = -g */
      if ( Parm->PrintLevel >= 1 ) printf ("RESTART CG\n") ;

      cg_precon ( gtemp, gtemp_P, _Com ) ;

      _gnorm    = cg_linf ( gtemp, n) ;
      _gnorm2_P = cg_dot ( gtemp, gtemp_P, n) ;
      
      cg_copy (x, xtemp, n) ;
      cg_copy (g, gtemp, n) ;

      for ( i = 0; i < n; i++ )
	d [i] = -gtemp_P[i] ;
	
      if ( cg_tol (_f, _gnorm, _StopRule, _tol) )
	{
	  *status_code = 0 ;
	  return 0 ;
	}

      _dphi0 = -_gnorm2_P ;
    }


  if ( !_Com->AWolfe )
    if ( fabs (_f-_Com->f0) < Parm->AWolfeFac*_Ck ) 
      _Com->AWolfe = TRUE ;

  if ( Parm->PrintLevel >= 1 )
#ifdef MPISRC
    ROOTONLY
#endif
      printf ("\niter: %5i f = %14.6e gnorm = %14.6e AWolfe = %2i\n",
	      (int) iter, _f, _gnorm, _Com->AWolfe) ;

  if ( Parm->debug )
    if ( _f > _Com->f0 + Parm->debugtol*_Ck )
      {
	*status_code = 9 ;
	return 0 ;
      }
  
  if ( _dphi0 > ZERO )
    {
      *status_code = 5 ;
      return 0 ;
    }
  return 1;
}

void CG_Descent::cg_post_process(double* x, int iter, cg_stats* Stat, int status_code)
{
  int i;
  double t;

  int status = status_code;
  int n = _Com->n;
  double* xtemp = _Com->xtemp;
  double* gtemp = _Com->gtemp;
  double* g     = _Com->g;
  cg_parameter* Parm = _Com->Parm;

  if ( Stat != NULL )
    {
      Stat->f     = _f ;
      Stat->gnorm = _gnorm ;
      Stat->nfunc = _Com->nf ;
      Stat->ngrad = _Com->ng ;
      Stat->iter  = iter ;
    }
  
  if ( status > 2 )
    {
      _gnorm = ZERO ;
      for (i = 0; i < n; i++)
	{
	  x [i] = xtemp [i] ;
	  g [i] = gtemp [i] ;
	  t = fabs (g [i]) ;
	  _gnorm = MAX (_gnorm, t) ;
        }

#ifdef MPISRC
      double tmp = DBL_MIN;
      gdmax(&_gnorm, 1, &tmp);
#endif

      if ( Stat != NULL ) Stat->gnorm = _gnorm ;
    }

#ifdef MPISRC
  ROOTONLY
#endif
    if ( Parm->PrintFinal || Parm->PrintLevel >= 1 )
      {
	const char mess1 [] = "Possible causes of this error message:" ;
	const char mess2 [] = "   - your tolerance may be too strict: "
	  "grad_tol = " ;
	const char mess3 [] = "Line search fails" ;
	const char mess4 [] = "   - your gradient routine has an error" ;
	const char mess5 [] = "   - the parameter epsilon in cg_descent_c.parm "
	  "is too small" ;
	printf ("\nTermination status: %i\n", status) ;
	if ( status == -2 )
	  {
	    printf ("At iteration %10.0f function value became nan\n",
		    (double) iter) ;
	  }
	else if ( status == -1 )
	  {
	    printf ("Objective function value is nan at starting point\n") ;
	  }
	else if ( status == 0 )
	  {
	    printf ("Convergence tolerance for gradient satisfied\n") ;
	  }
	else if ( status == 1 )
	  {
	    printf ("Terminating since change in function value "
		    "<= feps*|f|\n") ;
	  }
	else if ( status == 2 )
	  {
	    printf ("Number of iterations exceed specified limit\n") ;
	    printf ("Iterations: %10.0f maxit: %10.0f\n",
		    (double) iter, (double) _maxit) ;
	    printf ("%s\n", mess1) ;
	    printf ("%s %e\n", mess2, _grad_tol) ;
	  }
	else if ( status == 3 )
	  {
	    printf ("Slope always negative in line search\n") ;
	    printf ("%s\n", mess1) ;
	    printf ("   - your cost function has an error\n") ;
	    printf ("%s\n", mess4) ;
	  }
	else if ( status == 4 )
	  {
	    printf ("Line search fails, too many secant steps\n") ;
	    printf ("%s\n", mess1) ;
	    printf ("%s %e\n", mess2, _grad_tol) ;
	  }
	else if ( status == 5 )
	  {
	    printf ("Search direction not a descent direction\n") ;
	  }
	else if ( status == 6 ) /* line search fails */
	  {
	    printf ("%s\n", mess3) ;
	    printf ("%s\n", mess1) ;
	    printf ("%s %e\n", mess2, _grad_tol) ;
	    printf ("%s\n", mess4) ;
	    printf ("%s\n", mess5) ;
	  }
	else if ( status == 7 ) /* line search fails */
	  {
	    printf ("%s\n", mess3) ;
	    printf ("%s\n", mess1) ;
	    printf ("%s %e\n", mess2, _grad_tol) ;
	  }
	else if ( status == 8 ) /* line search fails */
	  {
	    printf ("%s\n", mess3) ;
	    printf ("%s\n", mess1) ;
	    printf ("%s %e\n", mess2, _grad_tol) ;
	    printf ("%s\n", mess4) ;
	    printf ("%s\n", mess5) ;
	  }
	else if ( status == 9 )
	  {
	    printf ("Debugger is on, function value does not improve\n") ;
	    printf ("new value: %25.16e old value: %25.16e\n", _f, _Com->f0) ;
	  }
	else if ( status == 10 )
	  {
	    printf ("Insufficient memory\n") ;
	  }
	
	printf ("maximum norm for gradient: %13.6e\n", _gnorm) ;
	printf ("function value:            %13.6e\n\n", _f) ;
	printf ("cg  iterations:          %10.0f\n", (double) iter) ;
	printf ("function evaluations:    %10.0f\n", (double) _Com->nf) ;
	printf ("gradient evaluations:    %10.0f\n", (double) _Com->ng) ;
	printf ("===================================\n\n") ;
      }

  return;
}

void CG_Descent::cg_precon(double* g, double* g_P, cg_com* Com)
{
  int n = Com->n ;

  if(!Com->cg_preconU)
    {
      cg_copy(g_P, g, n);
      return;
    }

  Com->cg_preconU(Com->instance, g, g_P, n, 1);

  return;
}

void CG_Descent::cg_precon_inv(double* d, double* d_invP, cg_com* Com)
{
  int n = Com->n ;
  
  if(!Com->cg_preconU)
    {
      cg_copy(d_invP, d, n);
      return ;
    }

  Com->cg_preconU(Com->instance, d, d_invP, n, 0);

  return;
}


/* =========================================================================
   === cg_default ==========================================================
   =========================================================================
   Set default conjugate gradient parameter values. If the parameter argument
   of cg_descent is NULL, this routine is called by cg_descent automatically.
   If the user wishes to set parameter values, then the cg_parameter structure
   should be allocated in the main program. The user could call cg_default
   to initialize the structure, and then individual elements in the structure
   could be changed, before passing the structure to cg_descent.
   =========================================================================*/
void CG_Descent::cg_default(cg_parameter* Parm)
{
  /* T => print final function value
     F => no printout of final function value */
  Parm->PrintFinal = TRUE ;

  /* Level 0 = no printing, ... , Level 3 = maximum printing */
  Parm->PrintLevel = 0 ;

  /* T => print parameters values
     F => do not display parmeter values */
  Parm->PrintParms = FALSE ;

  /* T => use approximate Wolfe line search
     F => use ordinary Wolfe line search, switch to approximate Wolfe when
     |f_k+1-f_k| < AWolfeFac*C_k, C_k = average size of cost */
  Parm->AWolfe = FALSE ;
  Parm->AWolfeFac = 1.e-3 ;

  /* factor in [0, 1] used to compute average cost magnitude C_k as follows:
     Q_k = 1 + (Qdecay)Q_k-1, Q_0 = 0,  C_k = C_k-1 + (|f_k| - C_k-1)/Q_k */
  Parm->Qdecay = .7 ;

  /* Stop Rules:
     T => ||grad||_infty <= max(grad_tol, initial |grad|_infty*StopFact)
     F => ||grad||_infty <= grad_tol*(1 + |f_k|) */
  Parm->StopRule = TRUE ;
  Parm->StopFac = 0.e-12 ;

  /* T => estimated error in function value is eps*Ck,
     F => estimated error in function value is eps */
  Parm->PertRule = TRUE ;
  Parm->eps = 1.e-6 ;

  /* T => attempt quadratic interpolation in line search when
     |f_k+1 - f_k|/f_k <= QuadCutoff
     F => no quadratic interpolation step */
  Parm->QuadStep = TRUE ;
  Parm->QuadCutOff = 1.e-12 ;

  /* T => check that f_k+1 - f_k <= debugtol*C_k
     F => no checking of function values */
  Parm->debug = FALSE ;
  Parm->debugtol = 1.e-10 ;

  /* if step is nonzero, it is the initial step of the initial line search */
  Parm->step = ZERO ;

  /* abort cg after maxit_fac*n iterations */
  Parm->maxit_fac = INF ;

  /* maximum number of times the bracketing interval grows or shrinks
     in the line search is nexpand */
  Parm->nexpand = (int) 50 ;

  /* maximum number of secant iterations in line search is nsecant */
  Parm->nsecant = (int) 50 ;

  /* conjugate gradient method restarts after (n*restart_fac) iterations */
  Parm->restart_fac = ONE ;

  /* stop when -alpha*dphi0 (estimated change in function value) <= feps*|f|*/
  Parm->feps = ZERO ;

  /* after encountering nan, growth factor when searching for
     a bracketing interval */
  Parm->nan_rho = 1.3 ;

  /* Wolfe line search parameter, range [0, .5]
     phi (a) - phi (0) <= delta phi'(0) */
  Parm->delta = .1 ;

  /* Wolfe line search parameter, range [delta, 1]
     phi' (a) >= sigma phi' (0) */
  Parm->sigma = .9 ;

  /* decay factor for bracket interval width in line search, range (0, 1) */
  Parm->gamma = .66 ;

  /* growth factor in search for initial bracket interval */
  Parm->rho = 5. ;

  /* conjugate gradient parameter beta_k must be >= eta*||d_k||_2 */
  Parm->eta = .01 ;

  /* starting guess for line search =
     psi0 ||x_0||_infty over ||g_0||_infty if x_0 != 0
     psi0 |f(x_0)|/||g_0||_2               otherwise */
  Parm->psi0 = .01 ;      /* factor used in starting guess for iteration 1 */

  /* for a QuadStep, function evalutated at psi1*previous step */
  Parm->psi1 = .1 ;

  /* when starting a new cg iteration, our initial guess for the line
     search stepsize is psi2*previous step */
  Parm->psi2 = 2. ;
}

/* =========================================================================
   ==== cg_Wolfe ===========================================================
   =========================================================================
   Check whether the Wolfe or the approximate Wolfe conditions are satisfied
   ========================================================================= */
int CG_Descent::cg_Wolfe
(
 double   alpha, /* stepsize */
 double       f, /* function value associated with stepsize alpha */
 double    dphi, /* derivative value associated with stepsize alpha */
 cg_com    *Com  /* cg com */
 )
{
  if ( dphi >= Com->wolfe_lo )
    {

      /* test original Wolfe conditions */

      if ( f - Com->f0 <= alpha*Com->wolfe_hi )
        {
	  if ( Com->Parm->PrintLevel >= 2 )
            {
	      printf ("wolfe f: %14.6e f0: %14.6e dphi: %14.6e\n",
		      f, Com->f0, dphi) ;
            }
	  return (1) ;
        }
      /* test approximate Wolfe conditions */
      else if ( Com->AWolfe )
        {
	  if ( (f <= Com->fpert) && (dphi <= Com->awolfe_hi) )
            {
	      if ( Com->Parm->PrintLevel >= 2 )
                {
		  printf ("f: %14.6e fpert: %14.6e dphi: %14.6e awolf_hi: "
			  "%14.6e\n", f, Com->fpert, dphi, Com->awolfe_hi) ;
                }
	      return (1) ;
            }
        }
    }
  return (0) ;
}

/* =========================================================================
   ==== cg_f ===============================================================
   Evaluate the function
   =========================================================================*/
double CG_Descent::cg_f
(
 double   *x,
 cg_com *Com
 )
{
  double f ;
  f = Com->cg_value (Com->instance, x, Com->n) ;
  Com->nf++ ;
  return (f) ;
}

/* =========================================================================
   ==== cg_g ===============================================================
   Evaluate the gradient
   =========================================================================*/
void CG_Descent::cg_g(double* g, double* x, cg_com* Com)
{
  Com->cg_grad (Com->instance, g, x, Com->n) ;
  Com->ng++ ;
}

/* =========================================================================
   ==== cg_fg ==============================================================
   Evaluate the function and gradient
   =========================================================================*/
double CG_Descent::cg_fg(double* g, double* x, cg_com* Com)
{
  double f ;
  if (Com->cg_valgrad != NULL ) f = Com->cg_valgrad (Com->instance, g, x, Com->n) ;
  else
    {
      Com->cg_grad (Com->instance, g, x, Com->n) ;
      f = Com->cg_value (Com->instance, x, Com->n) ;
    }
  Com->nf++ ;
  Com->ng++ ;
  return (f) ;
}

/* =========================================================================
   ==== cg_tol =============================================================
   =========================================================================
   Check for convergence
   ========================================================================= */
int CG_Descent::cg_tol
(
 double         f, /* function value associated with stepsize */
 double     gnorm, /* gradient sup-norm */
 int     StopRule, /* T => |grad|_infty <=max (tol, |grad|_infty*StopFact)
		      F => |grad|_infty <= tol*(1+|f|)) */
 double       tol  /* tolerance */
 )
{
  if ( StopRule )
    {
      if ( gnorm <= tol ) return (1) ;
    }
  else if ( gnorm <= tol*(ONE + fabs (f)) ) return (1) ;
  return (0) ;
}

/* =========================================================================
   ==== cg_dot =============================================================
   =========================================================================
   Compute dot product of x and y, vectors of length n
   ========================================================================= */
double CG_Descent::cg_dot
(
 double *x, /* first vector */
 double *y, /* second vector */
 int     n /* length of vectors */
 )
{
  double t = 0.0;

#ifdef OMPSRC

  int nthreads = (int) get_param("NTHREADS");
#pragma omp parallel num_threads(nthreads) default(shared)
  {
#pragma omp for reduction(+:t)
    for(int i = 0; i < n; i++)
      t += x[i]*y[i];
  }
#else
  int i, n5 ;
  n5 = n % 5 ;
  for (i = 0; i < n5; i++) t += x [i]*y [i] ;
  for (; i < n; i += 5)
    {
      t += x [i]*y[i] + x [i+1]*y [i+1] + x [i+2]*y [i+2]
	+ x [i+3]*y [i+3] + x [i+4]*y [i+4] ;
    }
#endif

#ifdef MPISRC
  double tmp = 0.0;
  gdsum(&t,1,&tmp);
#endif

  return (t) ;
}

/* =========================================================================
   ==== cg_linf =============================================================
   =========================================================================
   Compute L_inf norm of x , vectors of length n
   ========================================================================= */

double CG_Descent::cg_linf
(
 double *x, /* input vector */
 int     n  /* length of x */
 )
{
  double vmax = DBL_MIN;

#ifdef OMPSRC
  int nthreads = (int) get_param("NTHREADS");
#pragma omp parallel num_threads(nthreads) default(shared)
  {
    double vmax_loc = DBL_MIN;
    double tp;
#pragma omp for nowait
    for(int i = 0; i < n; i++)
      {
	tp = fabs(x[i]);
	vmax_loc = vmax_loc > tp? vmax_loc : tp;	
      }
#pragma omp critical
    vmax = vmax > vmax_loc? vmax:vmax_loc;
  }  
#else
  int i;
  int n5 = n%5 ;
  
  for (i = 0; i < n5; i++) vmax = vmax > fabs(x[i])? vmax:fabs(x[i]) ;
  for (; i < n; i += 5)
    { 
      vmax   = vmax > fabs(x[i])?   vmax : fabs(x[i]);
      vmax   = vmax > fabs(x[i+1])? vmax : fabs(x[i+1]);
      vmax   = vmax > fabs(x[i+2])? vmax : fabs(x[i+2]);
      vmax   = vmax > fabs(x[i+3])? vmax : fabs(x[i+3]);
      vmax   = vmax > fabs(x[i+4])? vmax : fabs(x[i+4]);
    }
#endif
 
#ifdef MPISRC
  double tmp = 0.0 ;
  gdmax(&vmax, 1, &tmp ) ;
#endif

  return vmax ;
}

void CG_Descent::cg_dvsub(int n, double *x, double *y, double *z)
{
#ifdef OMPSRC
  int nthreads = (int) get_param("NTHREADS");
#pragma omp parallel num_threads(nthreads) default(shared)
  {
#pragma omp for nowait
    for(int i = 0; i < n; i++)
      z[i] = x[i] - y[i];
  }
#else
  for(int i = 0; i < n; i++)
    z[i] = x[i] - y[i];
#endif

  return ;
}


/* =========================================================================
   === cg_copy =============================================================
   =========================================================================
   Copy vector x into vector y
   ========================================================================= */
void CG_Descent::cg_copy
(
 double *y, /* output of copy */
 double *x, /* input of copy */
 int     n  /* length of vectors */
 )
{
#ifdef OMPSRC
  int nthreads = (int) get_param("NTHREADS");
#pragma omp parallel num_threads(nthreads) default(shared)
  {
#pragma omp for nowait
    for(int i = 0; i < n; i++)
      y[i] = x[i];
  }
#else
  int j, n10 ;
  n10 = n % 10 ;
  for (j = 0; j < n10; j++) y [j] = x [j] ;
  for (; j < n; j += 10)
    {
      y [j] = x [j] ;
      y [j+1] = x [j+1] ;
      y [j+2] = x [j+2] ;
      y [j+3] = x [j+3] ;
      y [j+4] = x [j+4] ;
      y [j+5] = x [j+5] ;
      y [j+6] = x [j+6] ;
      y [j+7] = x [j+7] ;
      y [j+8] = x [j+8] ;
      y [j+9] = x [j+9] ;
    }
#endif

  return;
}

/* =========================================================================
   ==== cg_step ============================================================
   =========================================================================
   Compute xtemp = x + alpha d
   ========================================================================= */
void CG_Descent::cg_step
(
 double *xtemp, /*output vector */
 double     *x, /* initial vector */
 double     *d, /* search direction */
 double  alpha, /* stepsize */
 int         n  /* length of the vectors */
 )
{
#ifdef OMPSRC
  int nthreads = (int) get_param("NTHREADS");
#pragma omp parallel num_threads(nthreads) default(shared)
  {
#pragma omp for nowait
    for(int i = 0; i < n; i++)
      xtemp[i] = x[i] + alpha*d[i];
  }
#else
  int n5, i ;
  n5 = n % 5 ;
  for (i = 0; i < n5; i++) xtemp [i] = x[i] + alpha*d[i] ;
  for (; i < n; i += 5)
    { 
      xtemp [i]   = x [i]   + alpha*d [i] ;
      xtemp [i+1] = x [i+1] + alpha*d [i+1] ;
      xtemp [i+2] = x [i+2] + alpha*d [i+2] ;
      xtemp [i+3] = x [i+3] + alpha*d [i+3] ;
      xtemp [i+4] = x [i+4] + alpha*d [i+4] ;
    }
#endif

  return;
}

/* =========================================================================
   ==== cg_line ============================================================
   =========================================================================
   Approximate Wolfe line search routine
   ========================================================================= */
int CG_Descent::cg_line
(
 double  dphi0, /* function derivative at starting point (alpha = 0) */
 cg_com   *Com  /* cg com structure */
 )
{
  int n, iter ;
  int i, nsecant, nshrink, ngrow, status ;
  double a, dphia, b, dphib, c, alpha, phi, dphi,
    a0, da0, b0, db0, width, fquad, rho, *x, *xtemp, *d, *gtemp ;
  cg_parameter *Parm ;

  Parm = Com->Parm ;
  if ( Parm->PrintLevel >= 1 ) printf ("Approximate Wolfe line search\n") ;
  alpha = Com->alpha ;
  phi = Com->f ;
  n = Com->n ;
  x = Com->x ;         /* current iterate */
  xtemp = Com->xtemp ; /* x + alpha*d */
  d = Com->d ;         /* current search direction */
  gtemp = Com->gtemp ; /* gradient at x + alpha*d */
  rho = Parm->rho ;
  cg_step (xtemp, x, d, alpha, n) ;
  cg_g (gtemp, xtemp, Com) ;
  dphi = cg_dot (gtemp, d, n) ;

  /* check if gradient is nan; if so, reduce stepsize */
  if ( dphi != dphi )
    {
      for (i = 0; i < Parm->nexpand; i++)
        {
	  alpha /= rho ;
	  cg_step (xtemp, x, d, alpha, n) ;
	  cg_g (gtemp, xtemp, Com) ;
	  dphi = cg_dot (gtemp, d, n) ;
	  if ( dphi == dphi ) break ;
        }
      if ( i == Parm->nexpand )
        {
	  status = -2 ;
	  goto Exit ;
        }
      rho = Parm->nan_rho ;
    }
 
  /*Find initial interval [a,b] such that dphia < 0, dphib >= 0,
    and phia <= phi0 + feps*fabs (phi0) */
 
  a = ZERO ;
  dphia = dphi0  ;
  ngrow = 0 ;
  nshrink = 0 ;
  while ( dphi < ZERO )
    {
      phi = cg_f (xtemp, Com) ;

      /* if quadstep in effect and quadratic conditions hold, check wolfe condition*/

      if ( Com->QuadOK )
        {
	  if ( ngrow == 0 ) fquad = MIN (phi, Com->f0) ;
	  if ( phi <= fquad )
            {
	      if ( Parm->PrintLevel >= 2 )
                {
#ifdef MPISRC
		  ROOTONLY
#endif 
		    printf ("alpha: %14.6e phi: %14.6e fquad: %14.6e\n",
			    alpha, phi, fquad) ;
                }
	      if ( cg_Wolfe (alpha, phi, dphi, Com) )
                {
		  status = 0 ;
		  goto Exit ;
                }
            }
        }
      if ( phi > Com->fpert )
        {
	  /* contraction phase, only break at termination or Secant step */
	  b = alpha ;
	  while ( TRUE )
            {
	      alpha = .5*(a+b) ;
	      nshrink++ ;
	      if ( nshrink > Parm->nexpand )
                {
		  status = 6 ;
		  goto Exit ;
                }
	      cg_step (xtemp, x, d, alpha, n) ;
	      cg_g (gtemp, xtemp, Com) ;
	      dphi = cg_dot (gtemp, d, n) ;
	      if ( dphi >= ZERO ) goto Secant ;
	      phi = cg_f (xtemp, Com) ;
	      if ( Parm->PrintLevel >= 2 )
                {
#ifdef MPISRC
		  ROOTONLY
#endif 
		    printf ("contract, a: %14.6e b: %14.6e alpha: %14.6e phi: "
			    "%14.6e dphi: %14.6e\n", a, b, alpha, phi, dphi) ;
                }
	      if ( Com->QuadOK && (phi <= fquad) )
                {
		  if ( cg_Wolfe (alpha, phi, dphi, Com) )
                    {
		      status = 0 ;
		      goto Exit ;
                    }
                }
	      if ( phi <= Com->fpert )
                {
		  a = alpha ;
		  dphia = dphi ;
                }
	      else
                {
		  b = alpha ;
                }
            }
        }

      /* expansion phase */

      a = alpha ;
      dphia = dphi ;
      ngrow++ ;
      if ( ngrow > Parm->nexpand )
        {
	  status = 3 ;
	  goto Exit ;
        }
      alpha = rho*alpha ;
      cg_step (xtemp, x, d, alpha, n) ;
      cg_g (gtemp, xtemp, Com) ;
      dphi = cg_dot (gtemp, d, n) ;
      if ( Parm->PrintLevel >= 2 )
        {
#ifdef MPISRC
	  ROOTONLY
#endif 
	    printf ("expand,   a: %14.6e alpha: %14.6e phi: "
		    "%14.6e dphi: %14.6e\n", a, alpha, phi, dphi) ;
        }
    }

 Secant:
  b = alpha ;
  dphib = dphi ;
  if ( Com->QuadOK )
    {
      phi = cg_f (xtemp, Com) ;
      if ( ngrow + nshrink == 0 ) fquad = MIN (phi, Com->f0) ;
      if ( phi <= fquad )
        {
	  if ( cg_Wolfe (alpha, phi, dphi, Com) )
            {
	      status = 0 ;
	      goto Exit ;
            }
        }
    }
  nsecant = Parm->nsecant ;
  for (iter = 1; iter <= nsecant; iter++)
    {
      if ( Parm->PrintLevel >= 2 )
        {
#ifdef MPISRC
	  ROOTONLY
#endif 
	    printf ("secant, a: %14.6e b: %14.6e da: %14.6e db: %14.6e\n",
		    a, b, dphia, dphib) ;
        }
      width = Parm->gamma*(b - a) ;
      if ( -dphia <= dphib ) alpha = a - (a-b)*(dphia/(dphia-dphib)) ;
      else                   alpha = b - (a-b)*(dphib/(dphia-dphib)) ;
      c = alpha ;
      a0 = a ;
      b0 = b ;
      da0 = dphia ;
      db0 = dphib ;
      status = cg_update (&a, &dphia, &b, &dphib, &alpha, &phi, &dphi, Com) ;
      if ( status >= 0 ) goto Exit ;
      else if ( status == -2 )
        {
	  if ( c == a )
            {
	      if ( dphi > da0 ) alpha = c - (c-a0)*(dphi/(dphi-da0)) ;
	      else              alpha = a ;
            }
	  else
            {
	      if ( dphi < db0 ) alpha = c - (c-b0)*(dphi/(dphi-db0)) ;
	      else              alpha = b ;
            }
	  if ( (alpha > a) && (alpha < b) )
            {
	      if ( Parm->PrintLevel >= 2 ) 
#ifdef MPISRC
		ROOTONLY
#endif 
		  printf ("2nd secant\n") ;
	      status = cg_update (&a, &dphia, &b, &dphib, &alpha, &phi,
				  &dphi, Com) ;
	      if ( status >= 0 ) goto Exit ;
            }
        }

      /* bisection iteration */

      if ( b-a >= width )
        {
	  alpha = .5*(b+a) ;
	  if ( Parm->PrintLevel >= 2 ) 
#ifdef MPISRC
	    ROOTONLY
#endif 
	      printf ("bisection\n") ;
	  status = cg_update (&a, &dphia, &b, &dphib, &alpha, &phi,
			      &dphi, Com) ;
	  if ( status >= 0 ) goto Exit ;
        }
      else if ( b <= a )
        {
	  status = 7 ;
	  goto Exit ;
        }
    }
  status = 4 ;

 Exit:
  Com->alpha = alpha ;
  Com->f = phi ;
  Com->df = dphi ;
  return (status) ;
}

/* =========================================================================
   ==== cg_lineW ===========================================================
   =========================================================================
   Ordinary Wolfe line search routine.
   This routine is identical to cg_line except that the function
   psi [a] = phi [a] - phi [0] - a*delta*dphi [0] is minimized instead of
   the function phi
   ========================================================================= */
int CG_Descent::cg_lineW
(
 double  dphi0, /* function derivative at starting point (alpha = 0) */
 cg_com   *Com  /* cg com structure */
 )
{
  int n, iter ;
  int i, nsecant, nshrink, ngrow, status ;
  double a, dpsia, b, dpsib, c, alpha, phi, dphi,
    a0, da0, b0, db0, width, fquad, rho, psi, dpsi,
    *x, *xtemp, *d, *gtemp ;
  cg_parameter *Parm ;

  Parm = Com->Parm ;
  if ( Parm->PrintLevel >= 1 ) printf ("Wolfe line search\n") ;
  alpha = Com->alpha ;
  phi = Com->f ;
  dphi = Com->df ;
  n = Com->n ;
  x = Com->x ;         /* current iterate */
  xtemp = Com->xtemp ; /* x + alpha*d */
  d = Com->d ;         /* current search direction */
  gtemp = Com->gtemp ; /* gradient at x + alpha*d */
  rho = Parm->rho ;
  cg_step (xtemp, x, d, alpha, n) ;
  cg_g (gtemp, xtemp, Com) ;
  dphi = cg_dot (gtemp, d, n) ;

  /* check if gradient is nan; if so, reduce stepsize */
  if ( dphi != dphi )
    {
      for (i = 0; i < Parm->nexpand; i++)
        {
	  alpha /= rho ;
	  cg_step (xtemp, x, d, alpha, n) ;
	  cg_g (gtemp, xtemp, Com) ;
	  dphi = cg_dot (gtemp, d, n) ;
	  if ( dphi == dphi ) break ;
        }
      if ( i == Parm->nexpand )
        {
	  status = -2 ;
	  goto Exit ;
        }
      rho = Parm->nan_rho ;
    }
  dpsi = dphi - Com->wolfe_hi ;
 
  /*Find initial interval [a,b] such that dphia < 0, dphib >= 0,
    and phia <= phi0 + feps*fabs (phi0) */
 
  a = ZERO ;
  dpsia = dphi0 - Com->wolfe_hi ;
  ngrow = 0 ;
  nshrink = 0 ;
  while ( dpsi < ZERO )
    {
      phi = cg_f (xtemp, Com) ;
      psi = phi - alpha*Com->wolfe_hi ;

      /* if quadstep in effect and quadratic conditions hold, check Wolfe condition*/

      if ( Com->QuadOK )
        {
	  if ( ngrow == 0 ) fquad = MIN (phi, Com->f0) ;
	  if ( phi <= fquad )
            {
	      if ( Parm->PrintLevel >= 2 )
                {
#ifdef MPISRC
		  ROOTONLY
#endif 
		    printf ("alpha: %14.6e phi: %14.6e fquad: %14.6e\n",
			    alpha, phi, fquad) ;
                }
	      if ( cg_Wolfe (alpha, phi, dphi, Com) )
                {
		  status = 0 ;
		  goto Exit ;
                }
            }
        }
      if ( psi <= Com->fpert )
        {
	  a = alpha ;
	  dpsia = dphi ;
        }
      else
        {
	  /* contraction phase, only break at termination or Secant step */
	  b = alpha ;
	  while ( TRUE )
            {
	      alpha = .5*(a+b) ;
	      nshrink++ ;
	      if ( nshrink > Parm->nexpand )
                {
		  status = 6 ;
		  goto Exit ;
                }
	      cg_step (xtemp, x, d, alpha, n) ;
	      cg_g (gtemp, xtemp, Com) ;
	      dphi = cg_dot (gtemp, d, n) ;
	      dpsi = dphi - Com->wolfe_hi ;
	      if ( dpsi >= ZERO ) goto Secant ;
	      phi = cg_f (xtemp, Com) ;
	      psi = phi - alpha*Com->wolfe_hi ;
	      if ( Parm->PrintLevel >= 2 )
                {
#ifdef MPISRC
		  ROOTONLY
#endif 
		    printf ("contract, a: %14.6e b: %14.6e alpha: %14.6e phi: "
			    "%14.6e dphi: %14.6e\n", a, b, alpha, phi, dphi) ;
                }
	      if ( Com->QuadOK && (phi <= fquad) )
                {
		  if ( cg_Wolfe (alpha, phi, dphi, Com) )
                    {
		      status = 0 ;
		      goto Exit ;
                    }
                }
	      if ( psi <= Com->fpert )
                {
		  a = alpha ;
		  dpsia = dpsi ;
                }
	      else
                {
		  b = alpha ;
                }
            }
        }

      /* expansion phase */

      ngrow++ ;
      if ( ngrow > Parm->nexpand )
        {
	  status = 3 ;
	  goto Exit ;
        }
      alpha *= rho ;
      cg_step (xtemp, x, d, alpha, n) ;
      cg_g (gtemp, xtemp, Com) ;
      dphi = cg_dot (gtemp, d, n) ;
      dpsi = dphi - Com->wolfe_hi ;
      if ( Parm->PrintLevel >= 2 )
        {
#ifdef MPISRC
	  ROOTONLY
#endif 
	    printf ("expand,   a: %14.6e alpha: %14.6e phi: "
		    "%14.6e dphi: %14.6e\n", a, alpha, phi, dphi) ;
        }
    }

 Secant:
  b = alpha ;
  dpsib = dpsi ;
  if ( Com->QuadOK )
    {
      phi = cg_f (xtemp, Com) ;
      if ( ngrow + nshrink == 0 ) fquad = MIN (phi, Com->f0) ;
      if ( phi <= fquad )
        {
	  if ( cg_Wolfe (alpha, phi, dphi, Com) )
            {
	      status = 0 ;
	      goto Exit ;
            }
        }
    }
  nsecant = Parm->nsecant ;
  for (iter = 1; iter <= nsecant; iter++)
    {
      if ( Parm->PrintLevel >= 2 )
        {
#ifdef MPISRC
	  ROOTONLY
#endif 
	    printf ("secant, a: %14.6e b: %14.6e da: %14.6e db: %14.6e\n",
		    a, b, dpsia, dpsib) ;
        }
      width = Parm->gamma*(b - a) ;
      if ( -dpsia <= dpsib ) alpha = a - (a-b)*(dpsia/(dpsia-dpsib)) ;
      else                   alpha = b - (a-b)*(dpsib/(dpsia-dpsib)) ;
      c = alpha ;
      a0 = a ;
      b0 = b ;
      da0 = dpsia ;
      db0 = dpsib ;
      status = cg_updateW (&a, &dpsia, &b, &dpsib, &alpha, &phi, &dphi,
			   &dpsi, Com) ;
      if ( status >= 0 ) goto Exit ;
      else if ( status == -2 )
        {
	  if ( c == a )
            {
	      if ( dpsi > da0 ) alpha = c - (c-a0)*(dpsi/(dpsi-da0)) ;
	      else              alpha = a ;
            }
	  else
            {
	      if ( dpsi < db0 ) alpha = c - (c-b0)*(dpsi/(dpsi-db0)) ;
	      else              alpha = b ;
            }
	  if ( (alpha > a) && (alpha < b) )
            {
	      if ( Parm->PrintLevel >= 2 )
 #ifdef MPISRC
		ROOTONLY
#endif 
		  printf ("2nd secant\n") ;
	      status = cg_updateW (&a, &dpsia, &b, &dpsib, &alpha, &phi,
				   &dphi, &dpsi, Com) ;
	      if ( status >= 0 ) goto Exit ;
            }
        }

      /* bisection iteration */

      if ( b-a >= width )
        {
	  alpha = .5*(b+a) ;
	  if ( Parm->PrintLevel >= 2 ) 
#ifdef MPISRC
	    ROOTONLY
#endif 
	      printf ("bisection\n") ;
	  status = cg_updateW (&a, &dpsia, &b, &dpsib, &alpha, &phi, &dphi,
			       &dpsi, Com) ;
	  if ( status >= 0 ) goto Exit ;
        }
      else if ( b <= a )
        {
	  status = 7 ;
	  goto Exit ;
        }
    }
  status = 4 ;

 Exit:
  Com->alpha = alpha ;
  Com->f = phi ;
  Com->df = dphi ;
  return (status) ;
}

/* =========================================================================
   ==== cg_update ==========================================================
   =========================================================================
   update returns: 8 if too many iterations
   0 if Wolfe condition is satisfied
   -1 if interval is updated and a search is done
   -2 if the interval updated successfully
   ========================================================================= */
int CG_Descent::cg_update
(
 double        *a , /* left side of bracketing interval */
 double    *dphia , /* derivative at a */
 double        *b , /* right side of bracketing interval */
 double    *dphib , /* derivative at b */
 double    *alpha , /* trial step (between a and b) */
 double      *phi , /* function value at alpha (returned) */
 double     *dphi , /* function derivative at alpha (returned) */
 cg_com      *Com   /* cg com structure */
 )
{
  int n ;
  int nshrink, status ;
  double *x, *xtemp, *d, *gtemp ;
  cg_parameter *Parm ;

  Parm = Com->Parm ;
  n = Com->n ;
  x = Com->x ;         /* current iterate */
  xtemp = Com->xtemp ; /* x + alpha*d */
  d = Com->d ;         /* current search direction */
  gtemp = Com->gtemp ; /* gradient at x + alpha*d */
  cg_step (xtemp, x, d, *alpha, n) ;
  *phi = cg_fg (gtemp, xtemp, Com) ;
  *dphi = cg_dot (gtemp, d, n) ;
  if ( Parm->PrintLevel >= 2 )
    {
#ifdef MPISRC
      ROOTONLY
#endif 
	printf ("update alpha: %14.6e phi: %14.6e dphi: %14.6e\n",
		*alpha, *phi, *dphi) ;
    }
  if ( cg_Wolfe (*alpha, *phi, *dphi, Com) )
    {
      status = 0 ;
      goto Exit2 ;
    }
  status = -2 ;
  if ( *dphi >= ZERO )
    {
      *b = *alpha ;
      *dphib = *dphi ;
      goto Exit2 ;
    }
  else
    {
      if ( *phi <= Com->fpert )
        {
	  *a = *alpha ;
	  *dphia = *dphi ;
	  goto Exit2 ;
        }
    }
  nshrink = 0 ;
  *b = *alpha ;
  while ( TRUE )
    {
      *alpha = .5*(*a + *b) ;
      nshrink++ ;
      if ( nshrink > Parm->nexpand )
        {
	  status = 8 ;
	  goto Exit2 ;
        }
      cg_step (xtemp, x, d, *alpha, n) ;
      *phi = cg_fg (gtemp, xtemp, Com) ;
      *dphi = cg_dot (gtemp, d, n) ;
      if ( Parm->PrintLevel >= 2 )
        {
#ifdef MPISRC
	  ROOTONLY
#endif 
	    printf ("contract, a: %14.6e alpha: %14.6e "
		    "phi: %14.6e dphi: %14.6e\n", *a, *alpha, *phi, *dphi) ;
        }
      if ( cg_Wolfe (*alpha, *phi, *dphi, Com) )
        {
	  status = 0 ;
	  goto Exit2 ;
        }
      if ( *dphi >= ZERO )
        {
	  *b = *alpha ;
	  *dphib = *dphi ;
	  goto Exit1 ;
        }
      if ( *phi <= Com->fpert )
        {
	  if ( Parm->PrintLevel >= 2 )
            {
#ifdef MPISRC
	      ROOTONLY
#endif 
		printf ("update a: %14.6e dphia: %14.6e\n", *alpha, *dphi) ;
            }
	  *a = *alpha ;
	  *dphia = *dphi ;
        }
      else *b = *alpha ;
    }
 Exit1:
  status = -1 ;
 Exit2:
  if ( Parm->PrintLevel >= 2 )
    {
#ifdef MPISRC
      ROOTONLY
#endif 
	printf ("UP a: %14.6e b: %14.6e da: %14.6e db: %14.6e status: %i\n",
		*a, *b, *dphia, *dphib, status) ;
    }
  return (status) ;
}

/* =========================================================================
   ==== cg_updateW =========================================================
   =========================================================================
   This routine is identical to cg_update except that the function
   psi [a] = phi [a] - phi [0] - a*delta*dphi [0] is minimized instead of
   the function phi. The return int has the following meaning:
   8 if too many iterations
   0 if Wolfe condition is satisfied
   -1 if interval is updated and a search is done
   -2 if the interval updated successfully
   ========================================================================= */
int CG_Descent::cg_updateW
(
 double        *a , /* left side of bracketing interval */
 double    *dpsia , /* derivative at a */
 double        *b , /* right side of bracketing interval */
 double    *dpsib , /* derivative at b */
 double    *alpha , /* trial step (between a and b) */
 double      *phi , /* function value at alpha (returned) */
 double     *dphi , /* derivative of phi at alpha (returned) */
 double     *dpsi , /* derivative of psi at alpha (returned) */
 cg_com      *Com   /* cg com structure */
 )
{
  int n ;
  int nshrink, status ;
  double psi, *x, *xtemp, *d, *gtemp ;
  cg_parameter *Parm ;

  Parm = Com->Parm ;
  n = Com->n ;
  x = Com->x ;         /* current iterate */
  xtemp = Com->xtemp ; /* x + alpha*d */
  d = Com->d ;         /* current search direction */
  gtemp = Com->gtemp ; /* gradient at x + alpha*d */
  cg_step (xtemp, x, d, *alpha, n) ;
  *phi = cg_fg (gtemp, xtemp, Com) ;
  psi = *phi - *alpha*Com->wolfe_hi ;
  *dphi = cg_dot (gtemp, d, n) ;
  *dpsi = *dphi - Com->wolfe_hi ;
  if ( Parm->PrintLevel >= 2 )
    {
#ifdef MPISRC
      ROOTONLY
#endif 
	printf ("update alpha: %14.6e psi: %14.6e dpsi: %14.6e\n",
	      *alpha, psi, *dpsi) ;
    }
  if ( cg_Wolfe (*alpha, *phi, *dphi, Com) )
    {
      status = 0 ;
      goto Exit2 ;
    }
  status = -2 ;
  if ( *dpsi >= ZERO )
    {
      *b = *alpha ;
      *dpsib = *dpsi ;
      goto Exit2 ;
    }
  else
    {
      if ( psi <= Com->fpert )
        {
	  *a = *alpha ;
	  *dpsia = *dpsi ;
	  goto Exit2 ;
        }
    }
  nshrink = 0 ;
  *b = *alpha ;
  while ( TRUE )
    {
      *alpha = .5*(*a + *b) ;
      nshrink++ ;
      if ( nshrink > Parm->nexpand )
        {
	  status = 8 ;
	  goto Exit2 ;
        }
      cg_step (xtemp, x, d, *alpha, n) ;
      *phi = cg_fg (gtemp, xtemp, Com) ;
      *dphi = cg_dot (gtemp, d, n) ;
      *dpsi = *dphi - Com->wolfe_hi ;
      psi = *phi - *alpha*Com->wolfe_hi ;
      if ( Parm->PrintLevel >= 2 )
        {
#ifdef MPISRC
	  ROOTONLY
#endif 
	    printf ("contract, a: %14.6e alpha: %14.6e "
		    "phi: %14.6e dphi: %14.6e\n", *a, *alpha, *phi, *dphi) ;
        }
      if ( cg_Wolfe (*alpha, *phi, *dphi, Com) )
        {
	  status = 0 ;
	  goto Exit2 ;
        }
      if ( *dpsi >= ZERO )
        {
	  *b = *alpha ;
	  *dpsib = *dpsi ;
	  goto Exit1 ;
        }
      if ( psi <= Com->fpert )
        {
	  if ( Parm->PrintLevel >= 2 )
            {
#ifdef MPISRC
	      ROOTONLY
#endif 
		printf ("update a: %14.6e dpsia: %14.6e\n", *alpha, *dpsi) ;
            }
	  *a = *alpha ;
	  *dpsia = *dpsi ;
        }
      else *b = *alpha ;
    }
 Exit1:
  status = -1 ;
 Exit2:
  if ( Parm->PrintLevel >= 2 )
    {
#ifdef MPISRC
      ROOTONLY
#endif 
	printf ("UP a: %14.6e b: %14.6e da: %14.6e db: %14.6e status: %i\n",
		*a, *b, *dpsia, *dpsib, status) ;
    }
  return (status) ;
}

/* =========================================================================
   ==== cg_printParms ======================================================
   =========================================================================
   Print the contents of the cg_parameter structure
   ========================================================================= */
void CG_Descent::cg_printParms
(
 cg_parameter  *Parm
 )
{
  printf ("PARAMETERS:\n") ;
  printf ("\n") ;
  printf ("Wolfe line search parameter ..................... delta: %e\n",
	  Parm->delta) ;
  printf ("Wolfe line search parameter ..................... sigma: %e\n",
	  Parm->sigma) ;
  printf ("decay factor for bracketing interval ............ gamma: %e\n",
	  Parm->gamma) ;
  printf ("growth factor for bracket interval ................ rho: %e\n",
	  Parm->rho) ;
  printf ("growth factor for bracket interval after nan .. nan_rho: %e\n",
	  Parm->nan_rho) ;
  printf ("truncation factor for cg beta ..................... eta: %e\n",
	  Parm->eta) ;
  printf ("perturbation parameter for function value ......... eps: %e\n",
	  Parm->eps) ;
  printf ("factor for computing average cost .............. Qdecay: %e\n",
	  Parm->Qdecay) ;
  printf ("relative change in cost to stop quadstep ... QuadCufOff: %e\n",
	  Parm->QuadCutOff) ;
  printf ("factor multiplying gradient in stop condition . StopFac: %e\n",
	  Parm->StopFac) ;
  printf ("cost change factor, approx Wolfe transition . AWolfeFac: %e\n",
	  Parm->AWolfeFac) ;
  printf ("restart cg every restart_fac*n iterations . restart_fac: %e\n",
	  Parm->restart_fac) ;
  printf ("stop when cost change <= feps*|f| ................. eps: %e\n",
	  Parm->eps) ;
  printf ("starting guess parameter in first iteration ...... psi0: %e\n",
	  Parm->psi0) ;
  printf ("starting step in first iteration if nonzero ...... step: %e\n",
	  Parm->step) ;
  printf ("factor multiply starting guess in quad step ...... psi1: %e\n",
	  Parm->psi1) ;
  printf ("initial guess factor for general iteration ....... psi2: %e\n",
	  Parm->psi2) ;
  printf ("max iterations is n*maxit_fac ............... maxit_fac: %e\n",
	  Parm->maxit_fac) ;
  printf ("max expansions in line search ................. nexpand: %i\n",
	  Parm->nexpand) ;
  printf ("max secant iterations in line search .......... nsecant: %i\n",
	  Parm->nsecant) ;
  printf ("print level (0 = none, 2 = maximum) ........ PrintLevel: %i\n",
	  Parm->PrintLevel) ;
  printf ("Logical parameters:\n") ;
  if ( Parm->PertRule )
    printf ("    Error estimate for function value is eps\n") ;
  else
    printf ("    Error estimate for function value is eps*Ck\n") ;
  if ( Parm->QuadStep )
    printf ("    Use quadratic interpolation step\n") ;
  else
    printf ("    No quadratic interpolation step\n") ;
  if ( Parm->PrintFinal )
    printf ("    Print final cost and statistics\n") ;
  else
    printf ("    Do not print final cost and statistics\n") ;
  if ( Parm->PrintParms )
    printf ("    Print the parameter structure\n") ;
  else
    printf ("    Do not print parameter structure\n") ;
  if ( Parm->AWolfe)
    printf ("    Approximate Wolfe line search\n") ;
  else
    printf ("    Wolfe line search") ;
  if ( Parm->AWolfeFac > 0. )
    printf (" ... switching to approximate Wolfe\n") ;
  else
    printf ("\n") ;
  if ( Parm->StopRule )
    printf ("    Stopping condition uses initial grad tolerance\n") ;
  else
    printf ("    Stopping condition weighted by absolute cost\n") ;
  if ( Parm->debug)
    printf ("    Check for decay of cost, debugger is on\n") ;
  else
    printf ("    Do not check for decay of cost, debugger is off\n") ;
}


