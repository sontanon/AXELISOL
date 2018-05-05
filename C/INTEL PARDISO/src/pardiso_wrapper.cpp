#include "tools.h"
#include "pardiso_param.h"
#include "pardiso.h"

// Define for matrix, vector checks.
#undef DEBUG

void pardiso_wrapper(const csr_matrix A,// Matrix system to solve: Au = f.
	double *u,			// Solution array.
	double *f,			// RHS array.
	double *r,			// Residual, r = f - Au, array.
	const double tol,		// Tolerance convergence.
	double *norm,			// Pointer to final norm.
	int *convergence,		// Pointer to convergence flag.
	const int infnorm,		// Select infnorm or twonorm.
	const int lr_use,		// Low Rank update.
	const int perm_use,		// Calculate or use permutation.
					// 0: Do not use or calculate.
					// 1: Use permutation in perm array.
					// 2: Calculate permutation onto perm array.
	const int precond_use)		// Use previously computed LU with CGS iteration.
					// 0: Do not use CGS preconditioner.
					// L: Stopping criterion of Krylov-Subspace iteration 10**(-L).
{
	// Auxiliary doubles for residual.
	double res, res0;

	// Modify parameters according to permutation use.
	iparm[5 - 1] = perm_use;

	// Modify parameters according to CGS preconditioner.
	if (precond_use)
	{
		/// LU preconditioned with CGS.
		iparm[4 - 1] = 10 * precond_use + 1;
	}

// Debugging and recheck procedures.
#ifdef DEBUG
	// Check matrix for errors.
	iparm[27 - 1] = 1;
#endif

	// Reordering and symbolic factorization.
	phase = 11;
	pardiso(pt, &maxfct, &mnum, &mtype, &phase, 
		&n, A.a, A.ia, A.ja, perm, &nrhs, 
		iparm, &msglvl, &ddum, &ddum, &error);

	if (error != 0) 
	{
		printf("ERROR during symbolic factorization: %d.\n", error);
		exit(1);
	}
	
#ifdef VERBOSE
	printf("PARDISO: Reordering completed.\n");
	printf("PARDISO: Number of nonzeros in factors = %d.\n", iparm[18 - 1]);
	printf("PARDISO: Number of factorization MFLOPS = %d.\n", iparm[19 - 1]);
#endif

	// Numerical factorization.
	phase = 22;
	pardiso(pt, &maxfct, &mnum, &mtype, &phase, 
		&n, A.a, A.ia, A.ja, perm, &nrhs, 
		iparm, &msglvl, &ddum, &ddum, &error);

	if (error != 0) 
	{
		printf("ERROR during numerical factorization: %d.\n", error);
		exit(2);
	}

#ifdef VERBOSE
	printf("PARDISO: Factorization completed.\n");
#endif

	// Back substitution and iterative refinement.
	phase = 33;
	pardiso(pt, &maxfct, &mnum, &mtype, &phase, 
		&n, A.a, A.ia, A.ja, perm, &nrhs, 
		iparm, &msglvl, f, u, &error);

	// Report CGS iterations.
#ifdef VERBOSE
	if (precond_use)
		printf("PARDISO CGS PRECONDITIONER: iparm(20) = %d.\n", iparm[20 - 1]);
#endif

	if (error != 0) 
	{
		printf("ERROR during solution: %d,\n", error);
		exit(3);
	}


	// Compute residual with MKL CSR MV.
	struct matrix_descr descrA;
	sparse_matrix_t csrA;
	// Create hanlde with matrix.
	mkl_sparse_d_create_csr(&csrA, SPARSE_INDEX_BASE_ONE, A.nrows, A.ncols, A.ia, A.ia + 1, A.ja, A.a);
	// Create matrix description.
	descrA.type = SPARSE_MATRIX_TYPE_GENERAL;
	// Analyze sparse matrix: choose proper kernels and workload.
	mkl_sparse_optimize(csrA);
	// Compute r = alpha * A * u + beta * r.
	mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, -1.0, csrA, descrA, u, 0.0, r);
	// Add RHS.
	cblas_daxpy(A.nrows, 1.0, f, 1, r, 1);
	// Release memory.
	mkl_sparse_destroy(csrA);

	// Calculate norms.
	if (infnorm) 
	{
		res = ABS(r[cblas_idamax(A.nrows, r, 1)]);
		res0 = ABS(f[cblas_idamax(A.nrows, f, 1)]);
	}
	else 
	{
		res = cblas_dnrm2(A.nrows, r, 1);
		res0 = cblas_dnrm2(A.nrows, f, 1);
	}
	// Relative residual.
	res0 = res / res0;

#ifdef VERBOSE
	printf("PARDISO: Relative residual = %e.\n", res0);
	printf("PARDISO: Absolute residual = %e.\n", res);
#endif

	// Check residual and output convergence type.
	if (res0 < tol) 
	{
		// Relative convergence.
		*norm = res;
		*convergence = 1;
#ifdef VERBOSE
		printf("PARDISO: Converged relatively.\n");
#endif
	}
	else 
	{
		// No convergence.
		*norm = res;
		*convergence = 0;
#ifdef VERBOSE
		printf("\nPARDISO: WARNING: Failed to converge!\n\n");
#endif
	}

	// Return.
	return;
}
