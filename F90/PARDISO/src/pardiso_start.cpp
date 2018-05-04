// PARDISO global parameters header.
#define PARDISO_MAIN_FILE
#include "pardiso_param.h"
#include "pardiso.h"

#undef DEBUG

// Initialize PARDISO parameters and memory.
extern "C" void pardiso_start_(const int NrInterior, const int NzInterior)
{
	// Real unsymmetric matrix.
	mtype = 11;
	// One RHS.
	nrhs = 1;
	// Matrix dimension.
	n = (NrInterior + 2) * (NzInterior + 2);
	// Direct solver.
	solver = 0;
	// Auxiliary integer.
	int k = 0;

	// Set everything to zero beforehand.
	for (k = 0; k < 64; k++)
	{
		iparm[k] = 0;
		dparm[k] = 0.0;
		pt[k] = 0;
	}

#ifdef VERBOSE
	printf("PARDISO: Intializing via pardisoinit.\n");
#endif
	// Default setup via pardisoinit.
	pardisoinit(pt, &mtype, &solver, iparm, dparm, &error);

	// Check for errors regarding authenitcation.
	if (error != 0)
	{
		switch (error)
		{
			case -10:
				printf("ERROR: No license file found.\n");
				exit(-10);
			case -11:
				printf("ERROR: License is expired.\n");
				exit(-11);
			case -12:
				printf("ERROR: Wrong username or hostname.\n");
				exit(-12);
			default:
				printf("ERROR: Unknown error %d.\n", error);
				exit(1);
		}
	}
#ifdef VERBOSE
	printf("PARDISO: License check was succesful.\n");
#endif

	// Get number of processors via OMP_NUM_THREADS.
	char *omp_string = getenv("OMP_NUM_THREADS");
	int np = 1;
	if (omp_string != NULL)
	{
		sscanf(omp_string, "%d", &np);
	}
	else 
	{
		printf("ERROR: OMP_NUM_THREADS is not set.\n");
		exit(1);
	}

#ifdef VERBOSE
	printf("PARDISO: OMP_NUM_THREADS = %d.\n", np);
#endif

	// Problem fine-tune parameters. Used for perm_use = precond_use = 0.
	// Index minus 1 is done to reference to FORTRAN form in manual.
	iparm[1 - 1] = 1;	// Do not use default parameters.
	iparm[2 - 1] = 2;	// Parallel fill-in reordering from METIS.
	iparm[3 - 1] = np;	// Set number of procs.
	iparm[4 - 1] = 0;	// No iterative-direct algorithm.
	iparm[5 - 1] = 0;	// No user fill-in reducing permutation.
	iparm[6 - 1] = 0;	// Do not write solution into RHS.
	iparm[7 - 1] = 0;	// Not in use.
	iparm[8 - 1] = 10;	// Max numbers of iterative refinement steps.
	iparm[9 - 1] = 0;	// Not in use.
	iparm[10 - 1] = 13;	// Perturb the pivot elements with 1E-13.
	iparm[11 - 1] = 1;	// Use nonsymmetric permutation and scaling MPS.
	iparm[12 - 1] = 0;	// No conjugate transposed/transpose solve.
	iparm[13 - 1] = 1;      // Maximum weighted matching algorithm is switched-on (default for non-symmetric).
	iparm[14 - 1] = 0;	// Output: Number of perturbed pivots.
	iparm[15 - 1] = 0;	// Not in use.
	iparm[16 - 1] = 0;	// Not in use.
	iparm[17 - 1] = 0;	// Not in use.
	iparm[18 - 1] = 0;	// No Output: Number of nonzeros in the factor LU.
	iparm[19 - 1] = 0;	// No Output: Mflops for LU factorization.
	iparm[20 - 1] = 0;      // Output: Numbers of CG Iterations.
	iparm[24 - 1] = 1;	// Parallel Numerical Factorization.
	iparm[25 - 1] = 1;	// Parallel Forward/Backward Solve.
	maxfct = 1;		// Maximum number of numerical factorizations.
	mnum = 1;		// Which factorization to use.
	msglvl = MESSAGE_LEVEL;	// Print statistical information in file.
	error = 0;		// Initialize error flag.

	// Allocate permutation vector.
	perm = (int *)malloc(n * sizeof(int));

#ifdef VERBOSE
	printf("PARDISO: Setup solver memory and parameters.\n");
#endif

	// Setup matrix-vector multiplication type.
	// Non-transposed, i.e. y = A*x.
	uplo[0] = 'N';

#ifdef DEBUG
	printf("\nPARDISO IPARM and DPARM parameters:\n");
	for (k = 0; k < 64; k++)
	{
		printf("iparm[%d] = %d\tdparm[%d] = %lf\n", k + 1, iparm[k], k + 1, dparm[k]);
	}
#endif

	return;
}
