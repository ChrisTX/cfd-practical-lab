#ifndef FIELDS_H_
#define FIELDS_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "real.h"

// Create a 1D field of length size
REAL*	create1Dfield(int size);

// Create a 2D field of size sizeX x sizeY
REAL**	create2Dfield(int sizeX, int sizeY);

// Create a 1D field of length size
REAL*	createVector(int len);

// Create a 2D field of size sizeX x sizeY
REAL**	createMatrix(int rows, int cols);

// Delete a 1D field
void 	destroy1Dfield(REAL* field);

// Delete a 2D field with sizeY rows
void 	destroy2Dfield(REAL** field, int sizeY);

// Delete a 1D field
void	destroyVector(REAL* vector);

// Delete a 2D field with rows many rows
void	destroyMatrix(REAL** matrix, int rows);

// Fill a 1D field of size size with value
void	fill1Dfield(REAL value, REAL* field, int size);

// Fill a 2D field of size sizeX x sizeY with value
void	fill2Dfield(REAL value, REAL** field, int sizeX, int sizeY);

// Component wise comparison of field1 and field2, both of size size
int		isEqual1Dfield(REAL* field1, REAL* field2, int size);

// Component wise comparison of field1 and field2, both of size size
int		isEqual2Dfield(REAL** field1, REAL** field2, int sizeX, int sizeY);

// Applies the function func to all elements in field, which is of size size
void	applyFunctionTo1Dfield(REAL (*func)(REAL), REAL* field, int size);

// Applies the function func to all elements in field, which is of size sizeX * sizeY
void	applyFunctionTo2Dfield(REAL (*func)(REAL), REAL** field, int sizeX, int sizeY);

// Print out the 1D field field of size size
void	print1Dfield(REAL* field, int size);

// Print out the 2D field field of size sizeX x sizeY
void	print2Dfield(REAL** field, int sizeX, int sizeY);

// Print out the 1D field vector of size len
void	printVector(REAL* vector, int len);

// Print out the 1D field field of size size
void	printMatrix(REAL** matrix, int rows, int cols);

// Write the 1D field field of size size in the file fileName
void	write1Dfield(const char* fileName, REAL* field, int size);

// Write the 2D field field of size sizeX x sizeY in the file fileName
void	write2Dfield(const char* fileName, REAL** field, int sizeX, int sizeY);

// Read a 1D field from the file fileName and write its size in size
REAL* 	read1Dfield(const char* fileName, int* size);

// Read a 2D field from the file fileName and write its sizes in sizeX/sizeY
REAL**	read2Dfield(const char* fileName, int* sizeX, int* sizeY);

// Write a 2D scalar field in a VTK at file fileName with description. The field is given by field and its sizes are sizeX and sizeY. The points are rectangularly spaced with distances dx / dy in x / y direction 
void	writeVTKfileFor2DscalarField(const char* fileName, const char* description, REAL** field, int sizeX, int sizeY, REAL dx, REAL dy);

// Write a 2D vector field in a VTK at file fileName with description. The field is given by two coordinates, fieldU and fieldV of size sizeX and sizeY. The points are rectangularly spaced with distances dx / dy in x / y direction 
void	writeVTKfileFor2DvectorField(const char* fileName, const char* description, REAL** fieldU, REAL** fieldV, int sizeX, int sizeY, REAL dx, REAL dy);


#endif /* FIELDS_H_ */
