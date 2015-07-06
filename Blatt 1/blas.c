#include "blas.h"

void	axpy(REAL alpha, REAL* x, REAL* y, int len) {
	for(int i = 0; i < len; ++i)
		y[i] += alpha * x[i];
}

REAL	dot(REAL* x, REAL* y, int len) {
	REAL res = 0;
	for(int i = 0; i < len; ++i)
		res += x[i] * y[i];
	return res;
}

REAL	nrm2(REAL* x, int len) {
	REAL res = 0;
	for(int i = 0; i < len; ++i)
		res += x[i] * x[i];
	return sqrt(res);
}

void	copy(REAL* x, REAL* y, int len) {
	for(int i = 0; i < len; ++i)
		y[i] = x[i];
}

void	scal(REAL alpha, REAL* x, int len) {
	for(int i = 0; i < len; ++i)
		x[i] *= alpha;
}

void	gemv(REAL alpha, REAL** A, REAL* x, REAL beta, REAL* y, int rows, int cols) {
	for(int i = 0; i < rows; ++i)
		y[i] = alpha * dot(A[i], x, cols) + beta * y[i];
}

void	scal2Dfield(REAL alpha, REAL** X, int sizeX, int sizeY) {
	for(int i = 0; i < sizeX; ++i)
		scal(alpha, X[i], sizeY);
}

void	axpy2Dfield(REAL alpha, REAL** X, REAL** Y, int sizeX, int sizeY) {
	for(int i = 0; i < sizeX; ++i)
		axpy(alpha, X[i], Y[i], sizeY);
}