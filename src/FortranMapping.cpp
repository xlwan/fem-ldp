/*
 * FortranMapping.cpp
 *
 * Author: Xiaoliang Wan
 */

#include "FortranMapping.h"

/******************* Translations for using Fortran version of Blas ***********************/

/// - BLAS level 1: Copy \a x to \a y
void dcopyC(const int& n, const double *x, const int& incx, double *y,
	    const int& incy) 
{
  F77NAME(dcopy)(n, x, incx, y, incy);
}

/// - BLAS level 1: y = alpha \a x plus \a y
void daxpyC(const int& n, const double& alpha, const double *x,
	    const int& incx, const double *y, const int& incy) 
{
  F77NAME(daxpy)(n, alpha, x, incx, y, incy);
}
/// - BLAS level 1: Swap \a x with  \a y
void dswapC(const int& n, double *x, const int& incx, double *y,
	    const int& incy) 
{
  F77NAME(dswap)(n, x, incx, y, incy);
}
/// - BLAS level 1: x = alpha \a x
void dscalC(const int& n, const double& alpha, double *x, const int& incx) 
{
  F77NAME(dscal)(n, alpha, x, incx);
}
/// - BLAS level 1: Plane rotation by c = cos(theta), s = sin(theta)
void drotC(const int& n, double *x, const int& incx, double *y,
	   const int& incy, const double& c, const double& s) 
{
  F77NAME(drot)(n, x, incx, y, incy, c, s);
}
/// - BLAS level 1: output =   \f$ x^T  y \f$
double ddotC(const int& n, const double *x, const int& incx, const double *y,
	     const int& incy) 
{
  return F77NAME(ddot)(n, x, incx, y, incy);
}
/// - BLAS level 1: output = \f$ ||x||_2 \f$
double dnrm2C(const int& n, const double *x, const int& incx) 
{
  return F77NAME(dnrm2)(n, x, incx);
}
/// - BLAS level 1: output = \f$ ||x||_1 \f$
double dasumC(const int& n, const double *x, const int& incx) 
{
  return F77NAME(dasum)(n, x, incx);
}

/// - BLAS level 1: output = 1st value where \f$ |x[i]| = max |x|_1 \f$
/// Note it is modified to return a value between (0,n-1) as per
/// the standard C convention
int idamaxC(const int& n, const double *x, const int& incx) 
{
  return F77NAME(idamax)(n, x, incx) - 1;
}
/// - BLAS level 2: Matrix vector multiply y = A \e x where A[m x n]
void dgemvC(const char& trans, const int& m, const int& n, const double& alpha,
	    const double* a, const int& lda, const double* x, const int& incx,
	    const double& beta, double* y, const int& incy) 
{
  F77NAME(dgemv)(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
}

void dspmvC(const char& trans, const int& n, const double& alpha,
	    const double* a, const double* x, const int& incx, const double& beta,
	    double* y, const int& incy) 
{
  F77NAME(dspmv)(trans, n, alpha, a, x, incx, beta, y, incy);
}

void dgemmC(const char& transa, const char& transb, const int& m, const int& n,
	    const int& k, const double& alpha, const double* a, const int& lda,
	    const double* b, const int& ldb, const double& beta, double* c,
	    const int& ldc) 
{
  F77NAME(dgemm)(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

void Cdgemm(const int M, const int N, const int K, const double a, double *A,
	    const int ldA, double * B, const int ldB, const double b, double *C,
	    const int ldC) 
{
  dgemmC('N', 'N', N, M, K, a, B, ldB, A, ldA, b, C, ldC);
}

/******************* Translations for using Fortran version of Lapack ***********************/

/// factor a real packed-symmetric matrix using Bunch-Kaufman pivoting.
void dsptrfC(const char& uplo, const int& n, double* ap, int *ipiv, int& info) 
{
  F77NAME(dsptrf)(uplo, n, ap, ipiv, info);
}

/// Solve a real  symmetric matrix problem using Bunch-Kaufman pivoting.
void dsptrsC(const char& uplo, const int& n, const int& nrhs, const double* ap,
	     const int *ipiv, double* b, const int& ldb, int& info) 
{
  F77NAME(dsptrs)(uplo, n, nrhs, ap, ipiv, b, ldb, info);
}

/// Cholesky factor a real Positive Definite packed-symmetric matrix.
void dpptrfC(const char& uplo, const int& n, double *ap, int& info) 
{
  F77NAME(dpptrf)(uplo, n, ap, info);
}

/// Solve a real Positive defiinte symmetric matrix problem using
/// Cholesky factorization.
void dpptrsC(const char& uplo, const int& n, const int& nrhs, const double *ap,
	     double *b, const int& ldb, int& info) 
{
  F77NAME(dpptrs)(uplo, n, nrhs, ap, b, ldb, info);
}

/// Cholesky factorize a real positive-definite banded-symmetric matrix
void dpbtrfC(const char& uplo, const int& n, const int& kd, double *ab,
	     const int& ldab, int& info) 
{
  F77NAME(dpbtrf)(uplo, n, kd, ab, ldab, info);
}

/// Solve a real, Positive definite banded symmetric matrix problem
/// using Cholesky factorization.
void dpbtrsC(const char& uplo, const int& n, const int& kd, const int& nrhs,
	     const double *ab, const int& ldab, double *b, const int& ldb, int& info) 
{
  F77NAME(dpbtrs)(uplo, n, kd, nrhs, ab, ldab, b, ldb, info);
}

/// General banded matrix LU factorisation
void dgbtrfC(const int& m, const int& n, const int& kl, const int& ku,
	     double* a, const int& lda, int* ipiv, int& info) 
{
  F77NAME(dgbtrf)(m, n, kl, ku, a, lda, ipiv, info);
}

/// Solve general banded matrix using LU factorisation
void dgbtrsC(const char& trans, const int& n, const int& kl, const int &ku,
	     const int& nrhs, const double* a, const int& lda, int* ipiv, double* b,
	     const int& ldb, int& info) 
{
  F77NAME(dgbtrs)(trans, n, kl, ku, nrhs, a, lda, ipiv, b, ldb, info);
}

/// General matrix LU factorisation
void dgetrfC(const int& m, const int& n, double *a, const int& lda, int *ipiv,
	     int& info) 
{
  F77NAME(dgetrf)(m, n, a, lda, ipiv, info);
}

/// General matrix LU backsolve
void dgetrsC(const char& trans, const int& n, const int& nrhs, const double* a,
	     const int& lda, int* ipiv, double* b, const int& ldb, int& info) 
{
  F77NAME(dgetrs)(trans, n, nrhs, a, lda, ipiv, b, ldb, info);
}

/// Generate matrix inverse
void dgetriC(const int& n, double *a, const int& lda, const int *ipiv,
	     double *wk, const int& lwk, int& info) 
{
  F77NAME(dgetri)(n, a, lda, ipiv, wk, lwk, info);
}

//  -- Find eigenvalues of symmetric tridiagonal matrix
void dsterfC(const int& n, double *d, double *e, int& info) {
	F77NAME(dsterf)(n, d, e, info);
}

// -- Solve general real matrix eigenproblem.
void dgeevC(const char& uplo, const char& lrev, const int& n, double* a,
	    const int& lda, double* wr, double* wi, double* rev, const int& ldr,
	    double* lev, const int& ldv, double* work, const int& lwork, int& info) 
{
  F77NAME(dgeev)(uplo, lrev, n, a, lda, wr, wi, rev, ldr, lev, ldv, work,
		 lwork, info);
}

// -- Solve packed-symmetric real matrix eigenproblem.
void dspevC(const char& jobz, const char& uplo, const int& n, double* ap,
	    double* w, double* z, const int& ldz, double* work, int& info) 
{
  F77NAME(dspev)(jobz, uplo, n, ap, w, z, ldz, work, info);
}

// -- Solve packed-banded real matrix eigenproblem.
void dsbevC(const char& jobz, const char& uplo, const int& kl, const int& ku,
	    double* ap, const int& lda, double* w, double* z, const int& ldz,
	    double* work, int& info) 
{
  F77NAME(dsbev)(jobz, uplo, kl, ku, ap, lda, w, z, ldz, work, info);
}
