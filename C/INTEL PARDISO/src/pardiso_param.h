// Verbose output in pardiso_wrapper.
#undef VERBOSE
// PARDISO message level.
#define MESSAGE_LEVEL 0
// Standard headers for allocation.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// OPENBLAS or MKL implementation.
#define MKL
#ifdef MKL
#include "mkl_cblas.h"
#include "mkl_spblas.h"
#else
#include "cblas.h"
#endif

#ifdef PARDISO_MAIN_FILE
// Solver type.
int solver;
// Matrix type.
int mtype;
// Number of RHS.
int nrhs;
// Matrix dimension.
int n;
// Internal solver memory pointer pt.
void *pt[64];
// PARDISO control parameters.
int iparm[64];
// Maximum number of matrices to keep in memory.
int maxfct;
// Selected matrix to solve.
int mnum;
// PARDISO phase.
int phase;
// PARDISO error flag.
int error;
// PARDISO message level;
int msglvl;
// Dumb auxiliary variables.
double ddum;
int idum;
// Permutation vector.
int *perm;
// Low-rank vector.
int *diff;
// Matrix-vector multiplication type.
char uplo[1];
#else
extern int solver;
extern int mtype;
extern int nrhs;
extern int n;
extern void *pt[64];
extern int iparm[64];
extern int maxfct;
extern int mnum;
extern int phase;
extern int error;
extern int msglvl;
extern double ddum;
extern int idum;
extern int *perm;
extern int *diff;
extern char uplo[1];
#endif
