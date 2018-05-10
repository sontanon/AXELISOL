// Global headers files.
#include "tools.h"

// Elliptic solver headers.
#include "general_elliptic_csr_gen.h"
#include "pardiso_wrapper.h"
#include "elliptic_tools.h"

// Use infinity norm in solver.
#define INFNORM 0

#undef DEBUG

// General elliptic equation, solves the linear equation:
//     2       2       2          
// (a d  +  b d  +  c d  +  d d  +  e d  +  s) u = f,
//     rr      rz      zz      r       z
//
// where u is the solution and a, b, c, d, e, f, s are all 
// functions of (r, z).
// 
#ifdef FORTRAN
extern "C" void general_elliptic_(double *u,// output solution.
	double *res,		// output residual. 
	const double *ell_a,	// input a coefficient.
	const double *ell_b,	// input b coefficient.
	const double *ell_c,	// input c coefficient.
	const double *ell_d,	// input d coefficient.
	const double *ell_e,	// input e coefficient.
	const double *ell_s,	// input s coefficient.
	const double *ell_f,	// input f coefficient.
	const double *p_uInf,	 // u value at infinity for Robin BC.
	const int *p_robin,	 // Robin BC type: 1, 2, 3.
	const int *p_r_sym,	 // R symmetry: 1(even), -1(odd).
	const int *p_z_sym,	 // Z symmetry: 1(even), -1(odd).
	const int *p_NrInterior, // Number of r interior points.
	const int *p_NzInterior, // Number of z interior points.
	const int *p_ghost_zones,// Number of ghost zones.
	const double *p_dr, 	 // Spatial step in r.
	const double *p_dz,	 // Spatial step in z.
	const int *p_norder,	 // Finite difference evolution: 2 or 4.
	const int *p_lr_use,	 // Use low rank update.
	const int *p_precond_use)// Calculate and/or use preconditioner.
{
	// Variables passed by reference.
	double uInf = *p_uInf;
	int robin = *p_robin;
	int r_sym = *p_r_sym;
	int z_sym = *p_z_sym;
	int NrInterior = *p_NrInterior;
	int NzInterior = *p_NzInterior;
	int ghost_zones = *p_ghost_zones;
	double dr = *p_dr;
	double dz = *p_dz;
	int norder = *p_norder;
	int lr_use = *p_lr_use;
	int precond_use = *p_precond_use;
#else
void general_elliptic(double *u,// output solution.
	double *res,		// output residual. 
	const double *ell_a,	// input a coefficient.
	const double *ell_b,	// input b coefficient.
	const double *ell_c,	// input c coefficient.
	const double *ell_d,	// input d coefficient.
	const double *ell_e,	// input e coefficient.
	const double *ell_s,	// input s coefficient.
	const double *ell_f,	// input f coefficient.
	const double uInf,	// u value at infinity for robin bc.
	const int robin,	// robin bc type: 1, 2, 3.
	const int r_sym,	// r symmetry: 1(even), -1(odd).
	const int z_sym,	// z symmetry: 1(even), -1(odd).
	const int NrInterior,	// number of r interior points.
	const int NzInterior,	// number of z interior points.
	const int ghost_zones,	// number of ghost zones.
	const double dr,	// spatial step in r.
	const double dz,	// spatial step in z.
	const int norder,	// finite difference evolution: 2 or 4.
	const int lr_use,	// use low rank update.
	const int precond_use) 	// calculate and/or use preconditioner.
{
#endif
	// Set original number of ghost zones.
	int ghost = ghost_zones;

	// The main point of this solver is that it works on a smaller grid
	// than that used on the rest of the program.
	// For a second and fourth order approximations, we use a grid of 
	// NrInterior * NzInterior interior points plus a boundary of one point 
	// all arround it, thus a grid of (NrInterior + 2) * (NzInterior + 2).
	// 
	// Therefore a reduction is necessary, eliminating the lower-left sides 
	// of the grid which are later filled trivially using symmetry conditions.
	//
	// The current value of ghost zones is stored in temporary variable.
	// Ghost zones will later be reset to this original value.
	int temp_ghost = ghost;

	// Size of reduced arrays.
	size_t g_size = (NrInterior + 2) * (NzInterior + 2) * sizeof(double);

	// Allocate reduced arrays.
	double *g_u = (double *)malloc(g_size);
	double *g_f = (double *)malloc(g_size);
	double *g_a = (double *)malloc(g_size);
	double *g_b = (double *)malloc(g_size);
	double *g_c = (double *)malloc(g_size);
	double *g_d = (double *)malloc(g_size);
	double *g_e = (double *)malloc(g_size);
	double *g_s = (double *)malloc(g_size);
	double *g_res = (double *)malloc(g_size);

	// Reduce arrays.
	ghost_reduce(u, g_u, NrInterior, NzInterior, ghost);
	ghost_reduce(res, g_res, NrInterior, NzInterior, ghost);
	ghost_reduce(ell_a, g_a, NrInterior, NzInterior, ghost);
	ghost_reduce(ell_b, g_b, NrInterior, NzInterior, ghost);
	ghost_reduce(ell_c, g_c, NrInterior, NzInterior, ghost);
	ghost_reduce(ell_d, g_d, NrInterior, NzInterior, ghost);
	ghost_reduce(ell_e, g_e, NrInterior, NzInterior, ghost);
	ghost_reduce(ell_s, g_s, NrInterior, NzInterior, ghost);
	ghost_reduce(ell_f, g_f, NrInterior, NzInterior, ghost);

	// Set new ghost.
	ghost = 1;

	// Set temporary total number of points.
	int NrTotal = NrInterior + 2;
	int NzTotal = NzInterior + 2;

	// Allocate and generate CSR matrix.
	csr_matrix A;
	int DIM0 = NrTotal * NzTotal;
	int nnz0 = nnz_general_elliptic(NrInterior, NzInterior, norder, robin);
	csr_allocate(&A, DIM0, DIM0, nnz0);

	// Fill CSR matrix.
	csr_gen_general_elliptic(A, NrInterior, NzInterior, norder, dr, dz, g_a, g_b, g_c, g_d, g_e, g_s, g_f, uInf, robin, r_sym, z_sym);
	printf("GENREAL ELLIPTIC: Generated CSR matrix with %d rows, %d columns and %d nnz.\n", A.nrows, A.ncols, A.nnz);

	// Elliptic solver return variables.
	double norm = 0.0;
	int convergence = 0;
	double tol = (norder == 4) ? dr * dr * dz * dz : dr * dz;

	// Call elliptic solver.
	pardiso_wrapper(A, g_u, g_f, g_res, tol, &norm, &convergence, INFNORM, lr_use, precond_use);

	// Check solver convergence.
	if (convergence == 1)
	{
		printf("GENERAL ELLIPTIC: Solver converged!\n");
	}
	else
	{
		printf("GENERAL ELIPTIC: WARNING possible no convergence: %d.!\n", convergence);
	}
	printf("GENERAL ELLIPTIC: ||r|| = %3.3E.\n", norm);

	// Reset ghost and total number of points.
	ghost = temp_ghost;
	NrTotal = NrInterior + ghost + 1;
	NzTotal = NzInterior + ghost + 1;

	// Transfer solution and residual to original arrays.
	ghost_fill(g_u, u, r_sym, z_sym, NrInterior, NzInterior, ghost);
	ghost_fill(g_res, res, r_sym, z_sym, NrInterior, NzInterior, ghost);

	// Clear memory.
	free(g_u);
	free(g_res);
	free(g_a);
	free(g_b);
	free(g_c);
	free(g_d);
	free(g_e);
	free(g_s);
	free(g_f);

	// Clear CSR matrix.
	csr_deallocate(&A);

    return;
}
