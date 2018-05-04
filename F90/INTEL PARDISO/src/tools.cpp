// Headers.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <omp.h>
#include <string.h>

// CSR matrix index base.
#define BASE 1

// Indexing macro: requires that NzTotal be defined in scope.
#define IDX(i, j) ((i) * NzTotal + (j))

// CSR matrix type.
typedef struct csr_matrices
{
	// Nonzeros array.
	double *a;
	// Row pointer.
	int *ia;
	// Column pointer.
	int *ja;
	// Number of rows.
	int nrows;
	// Number of columns.
	int ncols;
	// Number of nonzeros.
	int nnz;

} csr_matrix;

// Grid parameters type.
typedef struct grid_params
{
	// Number of interior points in r.
	int NrInterior;
	// Number of interior points in z.
	int NzInterior;
	// Number of ghost zones.
	int ghost;
	// Finite difference order.
	int order;
	// Spatial step in r.
	double dr;
	// Spatial step in z.
	double dz;
} grid_param;

// Write simple ASCII 2D file, FORTRAN version.
extern "C" void write_single_file_(const double *u, const char *fname, const int NrTotal, const int NzTotal)
{
	// Auxiliary integers.
	int i, j;

	// Open file.
	FILE *fp = fopen(fname, "w");

	// Loop over r, z and write values.
	for (i = 0; i < NrTotal; i++)
	{
		for (j = 0; j < NzTotal; j++)
		{
			fprintf(fp, (j < NzTotal - 1) ? "%9.18E\t" : "%9.18E\n", u[IDX(i, j)]);
		}
	}

	// Close file.
	fclose(fp);

	return;
}

// C version.
void write_single_file(const double *u, const char *fname, const int NrTotal, const int NzTotal)
{
	// Auxiliary integers.
	int i, j;

	// Open file.
	FILE *fp = fopen(fname, "w");

	// Loop over r, z and write values.
	for (i = 0; i < NrTotal; i++)
	{
		for (j = 0; j < NzTotal; j++)
		{
			fprintf(fp, (j < NzTotal - 1) ? "%9.18E\t" : "%9.18E\n", u[IDX(i, j)]);
		}
	}

	// Close file.
	fclose(fp);

	return;
}

// Create CSR matrix.
void csr_allocate(csr_matrix *A, const int nrows, const int ncols, const int nnz)
{
	// Set integer parameters.
	A->nrows = nrows;
	A->ncols = ncols;
	A->nnz = nnz;
	// Allocate pointers.
	A->a = (double *)malloc(sizeof(double) * nnz);
	A->ja = (int *)malloc(sizeof(int) * nnz);
	A->ia = (int *)malloc(sizeof(int) * (nrows + 1));

	return;
}

// Destroy CSR matrix.
void csr_deallocate(csr_matrix *A)
{
	// Set integers to zero.
	A->nrows = 0;
	A->ncols = 0;
	A->nnz = 0;
	// Deallocate memory.
	free(A->a);
	free(A->ia);
	free(A->ja);

	return;
}

// CSR matrix print.
void csr_print(csr_matrix A, const char *vA, const char *iA, const char *jA)
{
	// Open three files corresponding to each array.
	FILE *fvA = fopen(vA, "w");
	FILE *fiA = fopen(iA, "w");
	FILE *fjA = fopen(jA, "w");

	// Auxiliary integers.
	int k = 0;

	// Loop over number of nonzeros and write A.a and A.ja.
	for (k = 0; k < A.nnz; k++)
	{
		fprintf(fvA, "%9.18E\n", A.a[k]);
		fprintf(fjA, "%d\n", A.ja[k]);
	}

	// Loop over A.nrows + 1 and write A.ia.
	for (k = 0; k < A.nrows + 1; k++)
		fprintf(fiA, "%d\n", A.ia[k]);

	// Close files.
	fclose(fvA);
	fclose(fiA);
	fclose(fjA);

	return;
}
