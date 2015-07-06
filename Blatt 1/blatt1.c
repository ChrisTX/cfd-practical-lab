#include "real.h"
#include "blas.h"
#include "fields.h"

#include <assert.h>

void testA() {
	REAL** A = create2Dfield(10, 10);
	assert(A);

	for(int i = 0; i < 10; ++i) {
		for(int j = 0; j < 10; ++j) {
			A[i][j] = (i+1)*(j+1);
		}
	}

	puts("A:");
	printMatrix(A, 10, 10);

	REAL* x = create1Dfield(10);
	REAL* y = create1Dfield(10);
	assert(x && y);

	for(int i = 0; i < 10; ++i) {
		x[i] = 1.;
		y[i] = 0.;
	}

	REAL* z = create1Dfield(10);
	gemv(1., A, x, 0., z, 10, 10);

	axpy(5., x, y, 10);
	axpy(1., z, y, 10);
	scal(1.5, x, 10);

	REAL dotvalue = dot(y, x, 10);
	REAL AXnorm = nrm2(z, 10);
	printf("nrm2(A*x): %f\n(5*x + A*x) * (1.5*x): %f\n", AXnorm, dotvalue);

	destroy1Dfield(x);
	destroy1Dfield(y);
	destroy2Dfield(A, 10);
}

void testB() {
	REAL** A = create2Dfield(60, 40);
	assert(A);

	const REAL pi = 2. * asin((REAL)1);

	for(int i = 0; i < 40; ++i) {
		for(int j = 0; j < 60; ++j)
			A[i][j] = pi * 1./60. * (REAL)i;
	}

	writeVTKfileFor2DscalarField("gradientField.vtk", "Gradient Field", A, 60, 40, (REAL)0.1, (REAL)0.1);

	#define SinGen(A) \
		_Generic(A, \
			float**: sinf, \
			double**: sin, \
			long double**: sinl)

	applyFunctionTo2Dfield(SinGen(A), A, 60, 40);

	writeVTKfileFor2DscalarField("sinusField.vtk", "Sined Gradient Field", A, 60, 40, (REAL)0.1, (REAL)0.1);

	destroy2Dfield(A, 40);
}

void testC() {
	int U_sizeX, U_sizeY;
	REAL** U = read2Dfield("fieldU.dat", &U_sizeX, &U_sizeY);
	assert(U);

	int V_sizeX, V_sizeY;
	REAL** V = read2Dfield("fieldV.dat", &V_sizeX, &V_sizeY);
	assert(V);

	assert(U_sizeX == V_sizeX && U_sizeY == V_sizeY);

	writeVTKfileFor2DvectorField("vectorField.vtk", "Vector Field", U, V, U_sizeX, U_sizeY, (REAL)0.1, (REAL)0.1);

	destroy2Dfield(V, V_sizeY);
	destroy2Dfield(U, U_sizeY);
}

int main() {
	// a)
	testA();

	// b)
	testB();

	// c)
	testC();
}