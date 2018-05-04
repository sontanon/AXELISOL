// Architecture: UNIX based vs. Windows.
#undef WIN

// Standard headers.
#include <stdio.h>
#include <stdlib.h>
#ifdef WIN
#define _USE_MATH_DEFINES
#endif
#include <math.h>
#include <time.h>
#include <assert.h>
#include <omp.h>
#include <string.h>

// System headers.
#include <sys/types.h>
#include <sys/stat.h>
#ifdef WIN
#include <io.h>
#include <direct.h>
#else
#include <unistd.h>
#endif

// Indexing macro: requires that NzTotal be defined in scope.
#define IDX(i, j) ((i) * NzTotal + (j))

// MIN/MAX macros.
#define MIN(X, Y) ((X) < (Y)) ? (X) : (Y)
#define MAX(X, Y) ((X) > (Y)) ? (X) : (Y)

// ABS macro.
#define ABS(X) ((X) < 0) ? -(X) : (X)

// CSR matrix index base.
#define BASE 1

// CSR matrix type.
typedef struct csr_matrices
{
	double *a;
	int *ia;
	int *ja;
	int nrows;
	int ncols;
	int nnz;

} csr_matrix;

// Forward declarations.
// 
// Write simple ASCII 2D file.
void write_single_file(const double *u, const char *fname, const int NrTotal, const int NzTotal);
// Create CSR matrix.
void csr_allocate(csr_matrix *A, const int nrows, const int ncols, const int nnz);
// Destroy CSR matrx.
void csr_deallocate(csr_matrix *A);
// CSR matrix print.
void csr_print(csr_matrix A, const char *vA, const char *iA, const char *jA);
