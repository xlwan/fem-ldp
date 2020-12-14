/*
 * Common.cpp
 *
 * 
 * Author: Xiaoliang Wan
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "Common.h"

static Param* param_list = 0;

void read_param(const char* fn)
{
  char   buf[BUFSIZ];
  char   name[NAMEMAX];
  double value;
  int    n;
  

  FILE* fp = fopen(fn, "r");
  fgets(buf, BUFSIZ, fp);
  fgets(buf, BUFSIZ, fp);

  Param* item;

  if(sscanf(buf, "%d", &n) != 1)
    {
      printf("ReadParams: cannot read # of parameters.\n");
      exit(-1);
    }

  while(n--)
    {
      fgets(buf, BUFSIZ, fp);
      item = new Param();
      sscanf(buf, "%lf%s", &item->val, item->name);
      item->next = param_list;
      param_list = item;
    }

  return;
}

Param* read_param(FILE* fp)
{
  char   buf[BUFSIZ];
  char   name[NAMEMAX];
  double value;
  int    n;

  fgets(buf, BUFSIZ, fp);
  fgets(buf, BUFSIZ, fp);

  Param* item;

  Param* list = 0;

  if(sscanf(buf, "%d", &n) != 1)
    {
      printf("ReadParams: cannot read # of parameters.\n");
      exit(-1);
    }

  while(n--)
    {
      fgets(buf, BUFSIZ, fp);
      item = new Param();
      sscanf(buf, "%lf%s", &item->val, item->name);
      item->next = list;
      list = item;
    }

  return list;
}


double get_param(const char* name)
{
  for(Param* p = param_list; p; p = p->next)
    if(strcmp(name, p->name) == 0)
      return p->val;
  printf("Fail to find the desired parameter: %s\n", name);
  exit(-1);
}

double get_param(const char* name, Param* list)
{
  for(Param* p = list; p; p = p->next)
    if(strcmp(name, p->name) == 0)
      return p->val;
  printf("Fail to find the desired parameter: %s\n", name);
  exit(-1);
}

void set_param(const char* name, double val)
{
  for(Param* p = param_list; p; p = p->next)
    if(strcmp(name, p->name) == 0)
      {
        p->val = val;
        return;
      }
  printf("Fail to find the desired parameter: %s\n", name);
  exit(-1);
}

void cs2hc(int n_cs, double* x_cs, int n_hc, double* x_hc) 
{
  
  int nhalf = n_hc > n_cs? n_cs:n_hc;
  
  nhalf /= 2;
  
  for(int i = 0; i < n_hc; i++) x_hc[i] = 0.0;
  
  x_hc[0] = x_cs[0];
  
  for(int i = 1; i <= nhalf; i++)
    x_hc[i] = x_cs[i]/2.0;
  
  for(int i = 1; i < nhalf; i++)
    x_hc[n_hc-i] = -x_cs[n_cs-i]/2.0;
}


void hc2cs(int n_hc, double* x_hc, int n_cs, double* x_cs)
{
  
  int nhalf = n_hc > n_cs? n_cs:n_hc;
  
  nhalf /= 2;
  
  for(int i = 0; i < n_cs; i++) x_cs[i] = 0.0;

  x_cs[0] = x_hc[0];
  
  for(int i = 1; i <= nhalf; i++)
    x_cs[i] = x_hc[i]*2.0;
  
  for(int i = 1; i < nhalf; i++)
    x_cs[n_cs-i] = -x_hc[n_hc-i]*2.0;
  
}

double ddotw(int n, double* x, int ix, double* y, int iy, double* w)
{
  double val = 0.0;

  while(n--)
    {
      val += (*x)*(*y)*(*w);
      x   += ix;
      y   += iy;
      w   += 1;
    }

  return val;
}

void dump2file(int M, int N, double* data, char* filename)
{
  FILE* fp = fopen(filename, "w");

  for(int i = 0; i < M; i++)
    {
      for(int j = 0; j < N; j++)
	fprintf(fp, "%20.12e ", data[i*N+j]);
      fprintf(fp,"\n");
    }

  fclose(fp);
}

void dvadd(int n, double* x, int incx, double* y, int incy, double* z, int incz)
{
  while(n--)
    {
      *z = (*x) + (*y);
      x += incx;
      y += incy;
      z += incz;
    }
  return;
}

void dvmul(int n, double* x, int incx, double* y, int incy, double* z, int incz)
{
  while(n--)
    {
      *z = (*x) * (*y);
      x += incx;
      y += incy;
      z += incz;
    }
  return;
}


void dvvtvp(int n, double *x, int incx, double *y, int incy,
                   double *w, int incw, double *z, int incz)
{
  while( n-- ) {
    *z = (*x) * (*y) + (*w);
    x += incx;
    y += incy;
    w += incw;
    z += incz;
  }
  return;
}


void dzero(int n, double* x, int incx)
{
  register double zero = 0;

  if(incx == 1)
    memset(x, '\0', n*sizeof(double));
  else
    while(n--)
      {
        *x  = zero;
         x += incx;
      }

  return;
}

double dsum(int n, double* x, int incx)
{
   double val = 0.0;
   
   while(n--)
     {
        val += *x;
        x   += incx;
     }

   return val;
}

double dmax(int n, double* x, int incx)
{
  double vmax = DBL_MIN;

  while(n--)
    {
      vmax = vmax > (*x)? vmax:(*x);
      x += incx;
    }
  return vmax;
}

double dmin(int n, double* x, int incx)
{
  double vmin = DBL_MAX;

  while(n--)
    {
      vmin = vmin < (*x)? vmin:(*x);
      x += incx;
    }

  return vmin;
}

double dmax_abs(int n, double* x, int incx)
{
  double vmax = 0;
  double tp;

  while(n--)
    {
      tp = fabs((*x));
      vmax = vmax > tp? vmax:tp;
      x += incx;
    }
  
  return vmax;
}

double dmin_abs(int n, double* x, int incx)
{
  double vmin = DBL_MAX;
  double tp;

  while(n--)
    {
      tp = fabs((*x));
      vmin = vmin < tp? vmin:tp;
      x += incx;
    }

  return vmin;
}


void dprint(int m, int n, double* x)
{
   for(int i = 0; i < m; i++)
     {
       for(int j = 0; j < n; j++)
         printf("%14.6e ", x[i*n+j]);
       printf("\n");
     }
    return;
}


int iwhere(double* p, int lb, int hb, double target)
{
  int iL = lb;
  int iH = hb;
  int iM;

  while(iH -iL > 1)
    {
      iM = (iH+iL)/2;
      if(target < p[iM])
	iH = iM;
      else if(target > p[iM])
	iL = iM;
      else
	return iM;
    }

  if(target == p[iH])
    return iH;

  return iL;
}


int iwhere_descend(double* p, int lb, int hb, double target)
{
  int iL = lb;
  int iH = hb;
  int iM;

  while(iH -iL > 1)
    {
      iM = (iH+iL)/2;
      if(target < p[iM])
	iL = iM;
      else if(target > p[iM])
	iH = iM;
      else
	return iM;
    }

  if(target == p[iH])
    return iH;

  return iL;
}


double** dmatrix(int m, int n)
{
  double** tp = new double*[m];
  tp[0]       = new double[m*n];

  for(int i = 0; i < m; i++)
    tp[i] = tp[0] + i*n;

  return tp;
}

int** imatrix(int m, int n)
{
  int** tp = new int*[m];
  tp[0]    = new int[m*n];

  for(int i = 0; i < m; i++)
    tp[i] = tp[0] + i*n;

  return tp;
}


void delete_matrix(double** p)
{
  delete[] p[0];
  delete[] p;
}

void delete_matrix(int** p)
{
  delete[] p[0];
  delete[] p;
}

double fs1(double s)
{
  double alpha = (double)get_param("ALPHA");

  double a = cos(M_PI - alpha);
  double b = cos(M_PI/2.0 + alpha);

  return (a+(b-a)*s);
}

double fs2(double s)
{
  double alpha = (double)get_param("ALPHA");

  double a = sin(M_PI - alpha);
  double b = sin(M_PI/2.0 + alpha);

  return (a+(b-a)*s);
}

double f1(double t)
{
  double T = 1.0;

  double val = -1.0 + 2.0*t/T;

  return val;
}

double f2(double t)
{
  double T = 1.0;
  double val;

  double v = 0.8;

  if(t <= T/4.0)
    {
      val = v*t/(T/4.0);
    }
  else if(t > T/4.0 && t < T*3.0/4.0)
    {
      val = v;
    }
  else
    {
      val = v - v/(T/4.0)*(t-T*3.0/4.0);
    }

  return val;
}

void tridag_solver(double* a, double* b, double* c, double* r, double* u, int N)
{
  double* gam = new double[N];
  double bet;
  
  bet  = b[0];
  u[0] = r[0]/bet;

  for(int i = 1; i < N; i++)
    {
      gam[i] = c[i-1]/bet;
      bet    = b[i] - a[i]*gam[i];
      if(!bet)
	{
	  printf("Error: tridag_solver failed!\n");
	  exit(-1);
	}
      u[i] = (r[i] - a[i]*u[i-1])/bet;
    }

  for(int i = N-2; i>=0; i--)
    u[i] = u[i] - gam[i+1]*u[i+1];

  delete[] gam;
  return;
}

// Following is a quicksort routine.
void swap(double*,double*);
int partition(double**, int, int);
int partition(double*, int, int);
int partition_descend(double**, int, int);
int partition_descend(double*, int, int);

void swap(double *left, double *right){
  double temp;
  
  temp  = *left;
  *left  = *right;
  *right = temp;
}

int partition(double **a, int left, int right){
  double target = a[0][right];
  int  i   = left -1;
  int  j   = right;
  
  while (1)
    {                                       
      while ( i < j )
	{
	  i++;
	  if (a[0][i] >= target) break;
	}
      
      while ( j > i )
	{
	  j--;
	  if (a[0][j] <= target) break;
	}
      
      if (i >= j)
	break;
      
      swap(a[0]+i, a[0]+j);
      swap(a[1]+i, a[1]+j);
    }
  
  swap(a[0]+i, a[0]+right);     
  swap(a[1]+i, a[1]+right);
  
  return i;
}

int partition_descend(double **a, int left, int right){
  double target = a[0][right];
  int  i   = left -1;
  int  j   = right;
  
  while (1)
    {                                       
      while ( i < j )
	{
	  i++;
	  if (a[0][i] <= target) break;
	}
      
      while ( j > i )
	{
	  j--;
	  if (a[0][j] >= target) break;
	}
      
      if (i >= j)
	break;
      
      swap(a[0]+i, a[0]+j);
      swap(a[1]+i, a[1]+j);
    }
  
  swap(a[0]+i, a[0]+right);     
  swap(a[1]+i, a[1]+right);
  
  return i;
}




int partition(double *a, int left, int right){
  double target = a[right];
  int  i   = left -1;
  int  j   = right;
  
  while (1)
    {                                       
      while ( i < j )
	{
	  i++;
	  if (a[i] >= target) break;
	}
      
      while ( j > i )
	{
	  j--;
	  if (a[j] <= target) break;
	}
      
      if (i >= j)
	break;
      
      swap(a+i, a+j);
    }
  
  swap(a+i, a+right);     
  
  return i;
}


int partition_descend(double *a, int left, int right){
  double target = a[right];
  int  i   = left -1;
  int  j   = right;
  
  while (1)
    {                                       
      while ( i < j )
	{
	  i++;
	  if (a[i] <= target) break;
	}
      
      while ( j > i )
	{
	  j--;
	  if (a[j] >= target) break;
	}
      
      if (i >= j)
	break;
      
      swap(a+i, a+j);
    }
  
  swap(a+i, a+right);     
  
  return i;
}




void quicksort(double** a, int left, int right)
{
  if(left >= right)
    return;

  int split = partition(a, left, right);
  
  quicksort(a, left,    split-1);
  quicksort(a, split+1, right  );       

  return;
}

void quicksort(double* a, int left, int right)
{
  if(left >= right)
    return;

  int split = partition(a, left, right);
  
  quicksort(a, left,    split-1);
  quicksort(a, split+1, right  );       

  return;
}

void quicksort_descend(double** a, int left, int right)
{
  if(left >= right)
    return;

  int split = partition_descend(a, left, right);
  
  quicksort_descend(a, left,    split-1);
  quicksort_descend(a, split+1, right  );       

  return;
}

void quicksort_descend(double* a, int left, int right)
{
  if(left >= right)
    return;

  int split = partition_descend(a, left, right);
  
  quicksort_descend(a, left,    split-1);
  quicksort_descend(a, split+1, right  );       

  return;
}

double timer()
{
  return (double)clock()/(double)CLOCKS_PER_SEC;
}

void exit_code(int status)
{
  exit(status);
}


