// Headers.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <omp.h>
#include <string.h>

// System headers.
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

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

// I/O subroutines in C for FORTRAN for GNU compatibility.
// Print help message.
#ifdef FORTRAN
extern "C" void print_help_(void)
{
	// User input char.
	char opt;

	printf("ELLSOLVEF: WARNING! Usage is  $./ELLSOLVEF dirname solver norder NrInterior NzInterior dr dz nrobin\n");
	printf("           [solver] is the type of solver: flat or general.\n");
	printf("           [dirname] is a valid directory string name.\n");
	printf("           [norder] is an integer equal to 2 or 4 corresponding to the finite difference order.\n");
	printf("           [NrInterior] and [NzInterior] are integers equal to the number of interior points in r, z.\n");
	printf("           [dr] and [dz] are floating point doubles equal to the spatial step in r, z.\n");
	printf("           [nrobin] is an optional argument corresponding to Robin operator order: 1, 2, 3.\n");
	printf("Press (y/n) to procede with default arguments:\n");
	opt = getchar();
	getchar();
	if ((opt =='y') || (opt == 'Y'))
	{
		printf("ELLSOLVEF: User chose to proceed with default arguments.\n");
	}
	else
	{
		printf("ELLSOLVEF: User chose to abort.\n");
		exit(1);
	}

	return;
}
#else
void print_help(void)
{
	// User input char.
	char opt;

	printf("ELLSOLVEC: WARNING! Usage is  $./ELLSOLVEC dirname solver norder NrInterior NzInterior dr dz nrobin\n");
	printf("           [solver] is the type of solver: flat or general.\n");
	printf("           [dirname] is a valid directory string name.\n");
	printf("           [norder] is an integer equal to 2 or 4 corresponding to the finite difference order.\n");
	printf("           [NrInterior] and [NzInterior] are integers equal to the number of interior points in r, z.\n");
	printf("           [dr] and [dz] are floating point doubles equal to the spatial step in r, z.\n");
	printf("           [nrobin] is an optional argument corresponding to Robin operator order: 1, 2, 3.\n");
	printf("Press (y/n) to procede with default arguments:\n");
	opt = getchar();
	getchar();
	if ((opt =='y') || (opt == 'Y'))
	{
		printf("ELLSOLVEC: User chose to proceed with default arguments.\n");
	}
	else
	{
		printf("ELLSOLVEC: User chose to abort.\n");
		exit(1);
	}

	return;
}
#endif

// Make directory and CD to it.
#ifdef FORTRAN
extern "C" void make_directory_and_cd_(const char *dirname)
#else
void make_directory_and_cd(const char *dirname)
#endif
{
	// User input char.
	char opt;

	// First create directory.
	struct stat st = { 0 };
	printf("I/O: Trying to create directory %s...\n", dirname);
	if (stat(dirname, &st) == -1)
	{

		mkdir(dirname, 0755);
	}
	else 
	{
		printf("I/O: WARNING! Directory %s already exists.\n", dirname);
		printf("Press (y/n) to procede and possibly overwrite files:\n");
		opt = getchar();
		getchar();
		if ((opt =='y') || (opt == 'Y'))
		{
			printf("I/O: User chose to proceed.\n");
		}
		else
		{
			printf("I/O: User chose to abort.\n");
			exit(1);
		}
	}
	// Now CD to output directory.
	if (chdir(dirname) == -1)
	{
		printf("I/O: WARNING! Could not CD to %s directory.\n", dirname);
		printf("Press (y/n) to procede and write in current directory:\n");
		opt = getchar();
		getchar();
		if ((opt =='y') || (opt == 'Y'))
		{
			printf("I/O: User chose to proceed.\n");
		}
		else
		{
			printf("I/O: User chose to abort.\n");
			exit(1);
		}
	}

	return;
}

// Write simple ASCII 2D file.
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
#ifdef FORTRAN
extern "C" void write_single_file_(const double *u, const char *fname, const int *p_NrTotal, const int *p_NzTotal)
{
	// Variables passed by reference.
	int NrTotal = *p_NrTotal;
	int NzTotal = *p_NzTotal;

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
#endif

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
