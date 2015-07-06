#include "fields.h"

REAL*	create1Dfield(int size) {
	return (REAL*)malloc(sizeof(REAL) * size);
}

REAL**	create2Dfield(int sizeX, int sizeY) {
	REAL** myrows = (REAL**)malloc(sizeof(REAL*) * sizeY);
	if(myrows == NULL)
		return myrows;

	int i = 0;
	int success = 1;
	for(; i < sizeY; ++i) {
		myrows[i] = create1Dfield(sizeX);
		if(myrows[i] == NULL) {
			success = 0;
			break;
		}
	}

	if(!success) {
		for(; i < sizeY; ++i)
			destroy1Dfield(myrows[i]);
		free(myrows);
		myrows = NULL;
	}

	return myrows;
}

REAL*	createVector(int len) {
	return create1Dfield(len);
}

REAL**	createMatrix(int rows, int cols) {
	return create2Dfield(rows, cols);
}

void 	destroy1Dfield(REAL* field) {
	free(field);
}

void 	destroy2Dfield(REAL** field, int sizeY) {
	for(int i = 0; i < sizeY; ++i)
		destroy1Dfield(field[i]);
	free(field);
}

void	destroyVector(REAL* vector) {
	destroy1Dfield(vector);
}

void	destroyMatrix(REAL** matrix, int rows) {
	destroy2Dfield(matrix, rows);
}

void	fill1Dfield(REAL value, REAL* field, int size) {
	for(int i = 0; i < size; ++i)
		field[i] = value;
}

void	fill2Dfield(REAL value, REAL** field, int sizeX, int sizeY) {
	for(int i = 0; i < sizeY; ++i)
		fill1Dfield(value, field[i], sizeX);
}

int		isEqual1Dfield(REAL* field1, REAL* field2, int size) {
	for(int i = 0; i < size; ++i)
		if(field1[i] != field2[i])
			return 0;

	return 1;
}

int		isEqual2Dfield(REAL** field1, REAL** field2, int sizeX, int sizeY) {
	for(int i = 0; i < sizeY; ++i)
		if(!isEqual1Dfield(field1[i], field2[i], sizeX))
			return 0;

	return 1;
}

void	applyFunctionTo1Dfield(REAL (*func)(REAL), REAL* field, int size) {
	for(int i = 0; i < size; ++i)
		field[i] = func(field[i]);
}

void	applyFunctionTo2Dfield(REAL (*func)(REAL), REAL** field, int sizeX, int sizeY) {
	for(int i = 0; i < sizeY; ++i)
		applyFunctionTo1Dfield(func, field[i], sizeX);
}

void	print1Dfield(REAL* field, int size) {
	for(int i = 0; i < size; ++i)
		printf("%f\n", field[i]);
}

void	print2Dfield(REAL** field, int sizeX, int sizeY) {
	for(int i = 0; i < sizeY; ++i) {
		for(int j = 0; j < sizeX; ++j) {
			if(j > 0)
				putchar(' ');

			printf("%f", field[i][j]);
		}
		putchar('\n');
	}
}

void	printVector(REAL* vector, int len) {
	print1Dfield(vector, len);
}

void	printMatrix(REAL** matrix, int rows, int cols) {
	print2Dfield(matrix, rows, cols);
}

void	write1Dfield(const char* fileName, REAL* field, int size) {
	FILE* myfile = fopen(fileName, "wb");
	if(myfile == NULL)
		return;

	int mysize = fwrite(&size, sizeof(size), 1, myfile);
	if(mysize == 1)
		fwrite(field, sizeof(REAL), size, myfile);

	fclose(myfile);
}

void	write2Dfield(const char* fileName, REAL** field, int sizeX, int sizeY) {
	FILE* myfile = fopen(fileName, "wb");
	if(myfile == NULL)
		return;

	int mysize = fwrite(&sizeX, sizeof(sizeX), 1, myfile);
	if(mysize < 1)
		goto write2Dend;

	mysize = fwrite(&sizeY, sizeof(sizeY), 1, myfile);
	if(mysize < 1)
		goto write2Dend;

	for(int i = 0; i < sizeY; ++i) {
		if(fwrite(field[i], sizeof(REAL), sizeX, myfile) != sizeX)
			goto write2Dend;
	}
	
write2Dend:
	fclose(myfile);
}

REAL* 	read1Dfield(const char* fileName, int* size) {
	FILE* myfile = fopen(fileName, "rb");
	REAL* retval = NULL;

	if(myfile == NULL)
		return NULL;

	if(!fread(size, sizeof(*size), 1, myfile))
		goto read1Dend;

	retval = create1Dfield(*size);
	if(!retval)
		goto read1Dend;

	if(fread(retval, sizeof(REAL), *size, myfile) != *size) {
		destroy1Dfield(retval);
		retval = NULL;
	}

read1Dend:
	fclose(myfile);
	return retval;
}

REAL**	read2Dfield(const char* fileName, int* sizeX, int* sizeY) {
	FILE* myfile = fopen(fileName, "rb");
	REAL** retval = NULL;

	if(myfile == NULL)
		return NULL;

	if(!fread(sizeX, sizeof(*sizeX), 1, myfile))
		goto read2Dend;

	if(!fread(sizeY, sizeof(*sizeY), 1, myfile))
		goto read2Dend;

	retval = create2Dfield(*sizeX, *sizeY);
	if(!retval)
		goto read2Dend;

	for(int i = 0; i < *sizeY; ++i) {
		if(fread(retval[i], sizeof(REAL), *sizeX, myfile) != *sizeX)
			goto read2Dkill;
	}

	goto read2Dend;

read2Dkill:
	destroy2Dfield(retval, *sizeY);
	retval = NULL;

read2Dend:
	fclose(myfile);
	return retval;
}

void	writeVTKfileFor2DscalarField(const char* fileName, const char* description, REAL** field, int sizeX, int sizeY, REAL dx, REAL dy) {
	FILE* myfile = fopen(fileName, "wb");

	if(myfile == NULL)
		return;

	int i, j;

	if(fprintf(myfile, "# vtk DataFile Version 3.0\nScalar Field\nASCII\nDATASET RECTILINEAR_GRID\nDIMENSIONS %d %d 1\nX_COORDINATES %d double\n", sizeX, sizeY, sizeX) < 0)
		goto write2Dend;

	for(i = 0; i < sizeX; ++i) {
		if(i > 0)
			if(fputc(' ', myfile) == EOF)
				goto write2Dend;

		if(fprintf(myfile, "%f", dx * i) < 0)
			goto write2Dend;
	}

	if(fprintf(myfile, "\nY_COORDINATES %d double\n", sizeY) < 0)
		goto write2Dend;

	for(i = 0; i < sizeY; ++i) {
		if(i > 0)
			if(fputc(' ', myfile) == EOF)
				goto write2Dend;

		if(fprintf(myfile, "%f", dy * i) < 0)
			goto write2Dend;
	}

	if(fprintf(myfile, "\nZ_COORDINATES 1 double\n0.0\nPOINT_DATA %d\nSCALARS %s double 1\nLOOKUP_TABLE default", sizeX * sizeY, description) < 0)
		goto write2Dend;

	for(i = 0; i < sizeX; ++i) {
		for(j = 0; j < sizeY; ++j) {
			if(fprintf(myfile, "\n%f", field[j][i]) < 0)
				goto write2Dend;
		}
	}

write2Dend:
	fclose(myfile);
}

void	writeVTKfileFor2DvectorField(const char* fileName, const char* description, REAL** fieldU, REAL** fieldV, int sizeX, int sizeY, REAL dx, REAL dy) {
	FILE* myfile = fopen(fileName, "wb");

	if(myfile == NULL)
		return;

	int i, j;

	if(fprintf(myfile, "# vtk DataFile Version 3.0\nVector Field\nASCII\nDATASET RECTILINEAR_GRID\nDIMENSIONS %d %d 1\nX_COORDINATES %d double\n", sizeX, sizeY, sizeX) < 0)
		goto write2Dvecend;

	for(i = 0; i < sizeX; ++i) {
		if(i > 0)
			if(fputc(' ', myfile) == EOF)
				goto write2Dvecend;

		if(fprintf(myfile, "%f", dx * i) < 0)
			goto write2Dvecend;
	}

	if(fprintf(myfile, "\nY_COORDINATES %d double\n", sizeY) < 0)
		goto write2Dvecend;

	for(i = 0; i < sizeY; ++i) {
		if(i > 0)
			if(fputc(' ', myfile) == EOF)
				goto write2Dvecend;

		if(fprintf(myfile, "%f", dy * i) < 0)
			goto write2Dvecend;
	}

	if(fprintf(myfile, "\nZ_COORDINATES 1 double\n0.0\nPOINT_DATA %d\nVECTORS %s double 1\nLOOKUP_TABLE default", sizeX * sizeY, description) < 0)
		goto write2Dvecend;

	for(i = 0; i < sizeX; ++i) {
		for(j = 0; j < sizeY; ++j) {
			if(fprintf(myfile, "\n%f %f 0.0", fieldU[j][i], fieldV[j][i]) < 0)
				goto write2Dvecend;
		}
	}

write2Dvecend:
	fclose(myfile);
}