// PARDISO prototypes for shared library.
extern "C" void pardisoinit (void *, int *, int *, int *, double *, int *);
extern "C" void pardiso(void *, int *, int *, int *, int *,
		int *, double *, int *, int *, int *, int *, 
		int *, int *, double *, double *, int *, double *);
extern "C" void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
extern "C" void pardiso_chkvec     (int *, int *, double *, int *);
extern "C" void pardiso_printstats (int *, int *, double *, int *, int *, int *, double *, int *);
