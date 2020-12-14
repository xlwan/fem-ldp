/*
 * FortranMapping.h
 *
 * interface between Fortran and C++ from commonly used Blas and Lapack routines.
 * 
 * Author: Xiaoliang Wan
 */

#ifndef FORTRANMAPPING_H_
#define FORTRANMAPPING_H_

#define F77NAME(x) x##_

/******************* Translations for using Fortran version of Blas ***********************/
extern "C" {
/// -- BLAS Level 1:
void F77NAME(dcopy)(const int& n, const double *x, const int& incx, double *y,
		const int& incy);
void F77NAME(daxpy)(const int& n, const double& alpha, const double *x,
		const int& incx, const double *y, const int& incy);
void F77NAME(dswap)(const int& n, double *x, const int& incx, double *y,
		const int& incy);
void F77NAME(dscal)(const int& n, const double& alpha, double *x,
		const int& incx);
void F77NAME(drot)(const int& n, double *x, const int& incx, double *y,
		const int& incy, const double& c, const double& s);
double F77NAME(ddot)(const int& n, const double *x, const int& incx,
		const double *y, const int& incy);
double F77NAME(dnrm2)(const int& n, const double *x, const int& incx);
double F77NAME(dasum)(const int& n, const double *x, const int& incx);
int F77NAME(idamax)(const int& n, const double *x, const int& incx);

/// -- BLAS level 2
void F77NAME(dgemv)(const char& trans, const int& m, const int& n,
		const double& alpha, const double* a, const int& lda, const double* x,
		const int& incx, const double& beta, double* y, const int& incy);

void F77NAME(dspmv)(const char& trans, const int& n, const double& alpha,
		const double* a, const double* x, const int& incx, const double& beta,
		double* y, const int& incy);

/// -- BLAS level 3:
void F77NAME(dgemm)(const char& trans, const char& transb, const int& m1,
		const int& n, const int& k, const double& alpha, const double* a,
		const int& lda, const double* b, const int& ldb, const double& beta,
		double* c, const int& ldc);
}

#endif /* FORTRANMAPPING_H_ */

/// - BLAS level 1: Copy \a x to \a y
void dcopyC(const int& n, const double *x, const int& incx, double *y,
		const int& incy);

/// - BLAS level 1: y = alpha \a x plus \a y
void daxpyC(const int& n, const double& alpha, const double *x,
		const int& incx, const double *y, const int& incy);

/// - BLAS level 1: Swap \a x with  \a y
void dswapC(const int& n, double *x, const int& incx, double *y,
		const int& incy);

/// - BLAS level 1: x = alpha \a x
void dscalC(const int& n, const double& alpha, double *x, const int& incx);

/// - BLAS level 1: Plane rotation by c = cos(theta), s = sin(theta)
void drotC(const int& n, double *x, const int& incx, double *y,
		const int& incy, const double& c, const double& s);

/// - BLAS level 1: output =   \f$ x^T  y \f$
double ddotC(const int& n, const double *x, const int& incx, const double *y,
		const int& incy);

/// - BLAS level 1: output = \f$ ||x||_2 \f$
double dnrm2C(const int& n, const double *x, const int& incx);

/// - BLAS level 1: output = \f$ ||x||_1 \f$
double dasumC(const int& n, const double *x, const int& incx);

/// - BLAS level 1: output = 1st value where \f$ |x[i]| = max |x|_1 \f$
/// Note it is modified to return a value between (0,n-1) as per
/// the standard C convention
int idamaxC(const int& n, const double *x, const int& incx);

/// - BLAS level 2: Matrix vector multiply y = A \e x where A[m x n]
void dgemvC(const char& trans, const int& m, const int& n, const double& alpha,
		const double* a, const int& lda, const double* x, const int& incx,
		const double& beta, double* y, const int& incy);

void dspmvC(const char& trans, const int& n, const double& alpha,
		const double* a, const double* x, const int& incx, const double& beta,
		double* y, const int& incy);

/// - BLAS level 3: Matrix-matrix multiply C = A x B where A[m x n],
///   B[n x k], C[m x k]
void dgemmC(const char& transa, const char& transb, const int& m, const int& n,
		const int& k, const double& alpha, const double* a, const int& lda,
		const double* b, const int& ldb, const double& beta, double* c,
		const int& ldc);

// Wrapper to mutliply two (row major) matrices together C = a*A*B + b*C
void Cdgemm(const int M, const int N, const int K, const double a, double *A,
		const int ldA, double * B, const int ldB, const double b, double *C,
		const int ldC);

/******************* Translations for using Fortran version of Lapack ***********************/
extern "C" {
// Matrix factorisation and solves
void F77NAME(dsptrf)(const char& uplo, const int& n, double* ap, int *ipiv,
		int& info);
void
		F77NAME(dsptrs)(const char& uplo, const int& n, const int& nrhs,
				const double* ap, const int *ipiv, double* b, const int& ldb,
				int& info);
void F77NAME(dpptrf)(const char& uplo, const int& n, double* ap, int& info);
void F77NAME(dpptrs)(const char& uplo, const int& n, const int& nrhs,
		const double* ap, double* b, const int& ldb, int& info);
void F77NAME(dpbtrf)(const char& uplo, const int& n, const int& kd, double* ab,
		const int& ldab, int& info);
void F77NAME(dpbtrs)(const char& uplo, const int& n, const int& kd,
		const int& nrhs, const double* ab, const int& ldab, double* b,
		const int& ldb, int& info);
void F77NAME(dgbtrf)(const int& m, const int& n, const int& kl, const int& ku,
		double* a, const int& lda, int* ipiv, int& info);
void F77NAME(dgbtrs)(const char& trans, const int& n, const int& kl,
		const int &ku, const int& nrhs, const double* a, const int& lda,
		int* ipiv, double* b, const int& ldb, int& info);
void F77NAME(dgetrf)(const int& m, const int& n, double* a, const int& lda,
		int* ipiv, int& info);
void F77NAME(dgetrs)(const char& trans, const int& n, const int& nrhs,
		const double* a, const int& lda, int* ipiv, double* b, const int& ldb,
		int& info);
void F77NAME(dgetri)(const int& n, double *a, const int& lda, const int *ipiv,
		double *wk, const int& lwk, int& info);
void F77NAME(dsterf)(const int& n, double *d, double *e, int& info);
void F77NAME(dgeev)(const char& uplo, const char& lrev, const int& n,
		double* a, const int& lda, double* wr, double* wi, double* rev,
		const int& ldr, double* lev, const int& ldv, double* work,
		const int& lwork, int& info);

void F77NAME(dspev)(const char& jobz, const char& uplo, const int& n,
		double* ap, double* w, double* z, const int& ldz, double* work,
		int& info);
void F77NAME(dsbev)(const char& jobz, const char& uplo, const int& kl,
		const int& ku, double* ap, const int& lda, double* w, double* z,
		const int& ldz, double* work, int& info);
}

/// factor a real packed-symmetric matrix using Bunch-Kaufman pivoting.
void dsptrfC(const char& uplo, const int& n, double* ap, int *ipiv, int& info);

/// Solve a real  symmetric matrix problem using Bunch-Kaufman pivoting.
void dsptrsC(const char& uplo, const int& n, const int& nrhs, const double* ap,
		const int *ipiv, double* b, const int& ldb, int& info);

/// Cholesky factor a real Positive Definite packed-symmetric matrix.
void dpptrfC(const char& uplo, const int& n, double *ap, int& info);

/// Solve a real Positive defiinte symmetric matrix problem using
/// Cholesky factorization.
void dpptrsC(const char& uplo, const int& n, const int& nrhs, const double *ap,
		double *b, const int& ldb, int& info);

/// Cholesky factorize a real positive-definite banded-symmetric matrix
void dpbtrfC(const char& uplo, const int& n, const int& kd, double *ab,
		const int& ldab, int& info);

/// Solve a real, Positive definite banded symmetric matrix problem
/// using Cholesky factorization.
void
		dpbtrsC(const char& uplo, const int& n, const int& kd, const int& nrhs,
				const double *ab, const int& ldab, double *b, const int& ldb,
				int& info);

/// General banded matrix LU factorisation
void dgbtrfC(const int& m, const int& n, const int& kl, const int& ku,
		double* a, const int& lda, int* ipiv, int& info);

/// Solve general banded matrix using LU factorisation
void dgbtrsC(const char& trans, const int& n, const int& kl, const int &ku,
		const int& nrhs, const double* a, const int& lda, int* ipiv, double* b,
		const int& ldb, int& info);

/// General matrix LU factorisation
void dgetrfC(const int& m, const int& n, double *a, const int& lda, int *ipiv,
		int& info);

/// General matrix LU backsolve
void dgetrsC(const char& trans, const int& n, const int& nrhs, const double* a,
		const int& lda, int* ipiv, double* b, const int& ldb, int& info);

/// Generate matrix inverse
void dgetriC(const int& n, double *a, const int& lda, const int *ipiv,
		double *wk, const int& lwk, int& info);

//  -- Find eigenvalues of symmetric tridiagonal matrix
void dsterfC(const int& n, double *d, double *e, int& info);

// -- Solve general real matrix eigenproblem.
void dgeevC(const char& uplo, const char& lrev, const int& n, double* a,
		const int& lda, double* wr, double* wi, double* rev, const int& ldr,
		double* lev, const int& ldv, double* work, const int& lwork, int& info);

// -- Solve packed-symmetric real matrix eigenproblem.

void dspevC(const char& jobz, const char& uplo, const int& n, double* ap,
		double* w, double* z, const int& ldz, double* work, int& info);

// -- Solve packed-banded real matrix eigenproblem.

void dsbevC(const char& jobz, const char& uplo, const int& kl, const int& ku,
		double* ap, const int& lda, double* w, double* z, const int& ldz,
		double* work, int& info);
