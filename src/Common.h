/*
 * Common.h
 *
 * Here is a collection of some useful commonly used routines.
 * 
 * Author: xlwan
 */


#ifndef COMMON_H_
#define COMMON_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <time.h>

#define NAMEMAX 50

typedef struct param
{
  char name[NAMEMAX];
  double val;
  param* next;
}Param;


// read the parameters form the input file.
void    read_param(const char* fn);
Param*  read_param(FILE* fp);

// get the value of a certain parameter from the Param list.
double  get_param(const char* name);
double  get_param(const char* name, Param* list);

void set_param(const char* name, double val);

/****************************************************************
 *
 *            Some small gadgets commonly used
 *
 ****************************************************************/

// transfer the cosine and sine coefficients into the halfcomplex format for FFT.
void cs2hc(int n_cs, double* x_cs, int n_hc, double* x_hc);

// transfer the halfcomplex data to cosine and sine coefficients.
void hc2cs(int n_hc, double* x_hc, int n_cs, double* x_nc);

// weighted dot
double ddotw(int n, double* x, int ix, double* y, int iy, double* w);

// dump an array to a file.
void dump2file(int Nx, int Ny, double* data, char* filename);

// summation of two vectors.
void dvadd(int n, double* x, int incx, double* y, int incy, double* z, int incz);

// multiplification of two vectors.
void dvmul(int n, double* x, int incx, double* y, int incy, double* z, int incz);

// z = x*y+w
void dvvtvp(int n, double* x, int incx, double* y, int incy, double* w, int incw, double* z, int incz);
void dzero(int n, double* x, int incx);

// sum of certain components of x.
double dsum(int n, double* x, int incx);

double dmax(int n, double* x, int incx);
double dmin(int n, double* x, int incx);
double dmax_abs(int n, double* x, int incx);
double dmin_abs(int n, double* x, int incx);

// print an array on the standard output.
void dprint(int m, int n, double* x);

// fast search of an arbitrarily given value in the array p.
int iwhere(double* p, int lb, int hb, double target);
int iwhere_descend(double* p, int lb, int hb, double target);

// allocate a matrix
double** dmatrix(int m, int n);
int** imatrix(int m, int n);

// free a matrix
void delete_matrix(double** p);
void delete_matrix(int** p);

// user given function.
double f1(double t);
double fs1(double s);

// user given function.
double f2(double t);
double fs2(double s);
  
// a solver for tridiagonal system
void tridag_solver(double* a, double* b, double* c, double* r, double* u, int N);

// quicksort
void quicksort(double** a, int left, int right);
void quicksort(double* a,  int left, int right);
void quicksort_descend(double** a, int left, int right);
void quicksort_descend(double* a, int left, int right);


double timer();
void exit_code(int status);


#endif /* COMMON_H_ */

