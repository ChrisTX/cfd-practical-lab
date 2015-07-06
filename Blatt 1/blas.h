#ifndef BLAS_H_
#define BLAS_H_

#include <math.h>
#include "real.h"

// Calculates y = alpha * x + y, where x and y are len long.
void	axpy(REAL alpha, REAL* x, REAL* y, int len);

// Calculates <x, y>, where x and y are len long.
REAL	dot(REAL* x, REAL* y, int len);

// Calculates ||x||_2, where x is len long.
REAL	nrm2(REAL* x, int len);

// Copies x to y, where x and y are len long.
void	copy(REAL* x, REAL* y, int len);

// Sets x = alpha * x, where x is len long.
void	scal(REAL alpha, REAL* x, int len);

// Calculates x = Ax + beta * y, where A is rows x cols, x is rows long and y is cols long
void	gemv(REAL alpha, REAL** A, REAL* x, REAL beta, REAL* y, int rows, int cols);

// Scales a sizeX x sizeY 2D field X by the factor alpha 
void	scal2Dfield(REAL alpha, REAL** X, int sizeX, int sizeY);

// Calculates Y = alpha * X + Y, where X and Y are sizeX x sizeY 2D fields
void	axpy2Dfield(REAL alpha, REAL** X, REAL** Y, int sizeX, int sizeY);
#endif /* BLAS_H_ */
