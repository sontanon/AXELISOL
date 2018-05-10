// PARDISO parameters and prototypes.
#include "pardiso_param.h"
#include "pardiso.h"

// Clear PARDISO memory and structures.
#ifdef FORTRAN
extern "C" void pardiso_stop_(void)
#else
void pardiso_stop(void)
#endif
{
#ifdef VERBOSE
	printf("PARDISO: Clearing internal memory...\n");
#endif
	// Termination and release of memory.
	phase = -1;
	pardiso(pt, &maxfct, &mnum, &mtype, &phase, 
		&n, &ddum, &idum, &idum, &idum, &nrhs,  
		iparm, &msglvl, &ddum, &ddum, &error);
#ifdef VERBOSE
	printf("PARDISO: Internal memory clear.\n");
#endif

	// Delete permutation vector.
	free(perm);
#ifdef VERBOSE
	printf("PARDISO: All memory clear.\n");
#endif

	return;
}
