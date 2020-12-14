#include <limits.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <ctype.h>

#define INT_INF INT_MAX
#define INF DBL_MAX

#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef NULL
#define NULL 0
#endif

#define ZERO ((double) 0)
#define ONE  ((double) 1)
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))


/*============================================================================
   user controlled parameters for the conjugate gradient algorithm
               (default values in cg_default)                                 */
typedef struct cg_parameter_struct /* user controlled parameters */
{
   /* parameters values that the user may wish to modify */
/*----------------------------------------------------------------------------*/
    /* T => print final statistics
       F => no printout of statistics */
    int PrintFinal ;

    /* Level 0  = no printing), ... , Level 3 = maximum printing */
    int PrintLevel ;

    /* T => print parameters values
       F => do not display parmeter values */
    int PrintParms ;

    /* T => use approximate Wolfe line search
       F => use ordinary Wolfe line search, switch to approximate Wolfe when
                |f_k+1-f_k| < AWolfeFac*C_k, C_k = average size of cost  */
    int    AWolfe ;
    double AWolfeFac ;

    /* factor in [0, 1] used to compute average cost magnitude C_k as follows:
       Q_k = 1 + (Qdecay)Q_k-1, Q_0 = 0,  C_k = C_k-1 + (|f_k| - C_k-1)/Q_k */
    double Qdecay ;

    /* Stop Rules:
       T => ||proj_grad||_infty <= max(grad_tol,initial ||grad||_infty*StopFact)
       F => ||proj_grad||_infty <= grad_tol*(1 + |f_k|) */
    int    StopRule ;
    double StopFac ;

    /* T => estimated error in function value is eps*Ck,
       F => estimated error in function value is eps */
    int    PertRule ;
    double eps ;

    /* T => attempt quadratic interpolation in line search when
                |f_k+1 - f_k|/f_k <= QuadCutoff
       F => no quadratic interpolation step */
    int    QuadStep ;
    double QuadCutOff ;

    /* T => check that f_k+1 - f_k <= debugtol*C_k
       F => no checking of function values */
    int    debug ;
    double debugtol ;

    /* if step is nonzero, it is the initial step of the initial line search */
    double step ;

    /* abort cg after maxit_fac*n iterations */
    double maxit_fac ;

    /* maximum number of times the bracketing interval grows or shrinks
       in the line search is nexpand */
    int nexpand ;

   /* maximum number of secant iterations in line search is nsecant */
    int nsecant ;

    /* conjugate gradient method restarts after (n*restart_fac) iterations */
    double restart_fac ;

    /* stop when -alpha*dphi0 (estimated change in function value) <= feps*|f|*/
    double feps ;

    /* after encountering nan, growth factor when searching for
       a bracketing interval */
    double nan_rho ;

/*============================================================================
       technical parameters which the user probably should not touch          */
    double           delta ; /* Wolfe line search parameter */
    double           sigma ; /* Wolfe line search parameter */
    double           gamma ; /* decay factor for bracket interval width */
    double             rho ; /* growth factor when searching for initial
                                bracketing interval */
    double             eta ; /* lower bound for the conjugate gradient update
                                parameter beta_k is eta*||d||_2 */
    double            psi0 ; /* factor used in starting guess for iteration 1 */
    double            psi1 ; /* in performing a QuadStep, we evaluate the
                                function at psi1*previous step */
    double            psi2 ; /* when starting a new cg iteration, our initial
                                guess for the line search stepsize is
                                psi2*previous step */
} cg_parameter ;

typedef struct cg_stats_struct /* statistics returned to user */
{
    double               f ; /*function value at solution */
    double           gnorm ; /* max abs component of gradient */
    int               iter ; /* number of iterations */
    int              nfunc ; /* number of function evaluations */
    int              ngrad ; /* number of gradient evaluations */
} cg_stats ;

// This is the core of control for cg_descent algorithms.
typedef struct cg_com_struct /* common variables */
{
  /* parameters computed by the code */
  int            n ; /* problem dimension, saved for reference */
  int           nf ; /* number of function evaluations */
  int           ng ; /* number of gradient evaluations */
  int         QuadOK ; /* T (quadratic step successful) */
  double       alpha ; /* stepsize along search direction */
  double           f ; /* function value for step alpha */
  double          df ; /* function derivative for step alpha */
  double       fpert ; /* perturbation is eps*Ck if PertRule is T */
  double          f0 ; /* old function value */
  double          Ck ; /* average cost as given by the rule:
			  Qk = Qdecay*Qk + 1, Ck += (fabs (f) - Ck)/Qk */
  double    wolfe_hi ; /* upper bound for slope in Wolfe test */
  double    wolfe_lo ; /* lower bound for slope in Wolfe test */
  double   awolfe_hi ; /* upper bound for slope, approximate Wolfe test */
  int         AWolfe ; /* F (use Wolfe line search)
			  T (use approximate Wolfe line search)
			  do not change user's AWolfe, this value can be
			  changed based on AWolfeFac */
  double          *x ; /* current iterate */
  double      *xtemp ; /* x + alpha*d */
  double          *d ; /* current search direction */
  double          *g ; /* current gradient */
  double      *gtemp ; /* gradient at x + alpha*d */
   double     *d_invP ; /* inv of P times d: inv(P)*d */
  double        *g_P ; /* Preconditioner times g: P*g*/
  double    *gtemp_P ; /* Preconditioner times gradient at x + alpha*d*/
  void      *instance; /* socket for a connection with outside */    	
  double   (*cg_value) (void *, double *, int) ;                /* f = cg_value (x, n) */
  void      (*cg_grad) (void *, double *, double *, int) ;      /* cg_grad (g, x, n) */
  double (*cg_valgrad) (void *, double *, double *, int) ;      /* f = cg_valgrad (g,x,n)*/
  void   (*cg_preconU) (void *, double *, double *, int, int) ; /* preconditioner routine by user */

  cg_parameter *Parm ; /* user parameters */
} cg_com ;

class CG_Descent{
 private:
 public:
  double _f;
  double _gnorm;
  double _gnorm2;

  CG_Descent();
  CG_Descent(int n);

  // allocate memory
  void cg_locate_mem(int n);

  // set parameter used by CG_Descent library.
  int cg_set_parameter(cg_parameter* UParm, void* instance, 
		       double (*value)  (void*, double*, int), 
		       void   (*grad)   (void*, double*, double*, int), 
		       double (*valgrad)(void*, double*, double*, int),
		       void   (*preconU)(void*, double*, double*, int, int));

  // initialize the iteration.
  int cg_initial_work(double* x, double grad_tol, int* maxit, int* status_code);
  
  // compute x_k+1 = x_k + alpha_k d_k
  int cg_forward_one(double* x, int* status_code);
  
  // compute the new direction: d_k+1 = g_k + beta_k d_k
  int cg_update_direction(double* x, int iter, int* status_code);
  
  // post processing.
  void cg_post_process(double* x, int iter, cg_stats* Stat, int status_code);

  /* set default parameter values */
  void cg_default(cg_parameter* Parm);

 protected:
  int    _maxit;
  int    _nrestart;
  int    _StopRule;
  double _dphi;
  double _dphi0;
  double _delta2;
  double _eta_sq;
  double _alpha;
  double _Ck;
  double _Qk;
  double _tol;
  double _grad_tol;
  double _gnorm_P;
  double _gnorm2_P;

  double* _work;
  cg_com* _Com;

  // g_P is equal to precoditioner times g.
  void cg_precon(double* g, double* g_P, cg_com* Com);

  // d_invP is equal to inverse of preconditioner times d.
  void cg_precon_inv(double* d, double* d_invP, cg_com* Com);
  
  int cg_Wolfe(double alpha, double f, double dphi, cg_com* Com);
  /*
    alpha: stepsize 
        f: function value associated with stepsize alpha 
     dphi: derivative value associated with stepsize alpha 
      Com: cg com 
   */

  // compute function
  double cg_f(double* x, cg_com* Com);

  // compute gradient
  void cg_g(double* g, double* x, cg_com* Com);

  // compute function and gradient
  double cg_fg(double* g, double* x, cg_com* Com);
  
  // check accuracy
  int cg_tol(double f, double gnorm, int StopRule, double tol);
  /*
           f: function value associated with stepsize 
       gnorm: gradient sup-norm 
    StopRule: T => |grad|_infty <=max (tol, |grad|_infty*StopFact)            
              F => |grad|_infty <= tol*(1+|f|)) 
         tol: tolerance 
  */

  // inner product of two vectors.
  double cg_dot(double* x, double* y, int  n);

  // compute the L_inf norm of a vector
  double cg_linf(double* x, int n);

  // vector operation: z = x - y
  void cg_dvsub(int n, double *x, double *y, double *z);

  // copy vector x to y.
  void cg_copy(double* y, double* x, int n); 

  // update one step of x
  void cg_step(double* xtemp, double* x, double* d, double alpha, int n);
  /*
    xtemp: output vector 
        x: initial vector 
        d: search direction 
    alpha: stepsize 
        n: length of the vectors 
   */
  
  // line search with approximate Wolfe conditions.
  int cg_line(double dphi0, cg_com* Com);
  /*
    dphi0: function derivative at starting point (alpha = 0) 
      Com: cg com structure 
  */

  // line search with Wolfe conditions.
  int cg_lineW(double dphi0, cg_com* Com);
  /*
    dphi0: function derivative at starting point (alpha = 0) 
      Com: cg com structure 
  */

  int cg_update(double* a, double* dphia, double* b, double* dphib, double* alpha, double* phi, double* dphi, cg_com* Com);
  /*
         a: left side of bracketing interval 
     dphia: derivative at a 
         b: right side of bracketing interval 
     dphib: derivative at b 
     alpha: trial step (between a and b) 
       phi: function value at alpha (returned) 
      dphi: function derivative at alpha (returned) 
       Com: cg com structure 
  */

  int cg_updateW(double* a, double* dpsia, double* b, double* dpsib, double* alpha, double* phi, double* dphi, double* dpsi, cg_com* Com);
  /*
       a: left side of bracketing interval 
   dpsia: derivative at a 
       b: right side of bracketing interval 
   dpsib: derivative at b 
   alpha: trial step (between a and b) 
     phi: function value at alpha (returned) 
    dphi: derivative of phi at alpha (returned) 
    dpsi: derivative of psi at alpha (returned) 
     Com: cg com structure 
   */

  void cg_printParms(cg_parameter* Parm);
};

