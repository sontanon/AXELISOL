void pardiso_wrapper(const csr_matrix A,// Matrix system to solve: Au = f.
	double *u,			// Solution array.
	double *f,			// RHS array.
	double *r,			// Residual, r = f - Au, array.
	const double tol,		// Tolerance convergence.
	double *norm,			// Pointer to final norm.
	int *convergence,		// Pointer to convergence flag.
	const int infnorm,		// Select infnorm or twonorm.
	const int precond_use);		// Use previously computed LU with CGS iteration.
					// 0: Do not use CGS preconditioner.
					// L: Stopping criterion of Krylov-Subspace iteration 10**(-L).
