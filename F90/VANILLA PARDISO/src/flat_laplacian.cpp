// Global header files.
#include "tools.h"

// Elliptic solver headers.
#include "flat_laplacian_csr_gen.h"
#include "pardiso_wrapper.h"
#include "elliptic_tools.h"

// Use infinity norm in solver: 0(twonorm), 1(infnorm).
#define INFNORM 0

#undef DEBUG

//  Flat Laplacian, solves the linear equation:
//    __2
//  ( \/     + s(r, z) ) u(r, z) = f(r, z).
//  	flat
//
//  Where the Laplacian is the flat Laplacian in cylindrical
//  coordinates:
//
//   __2       2      2
//   \/     = d   +  d    + (1/r) d  .
//     flat    rr     zz           r
//
//  And is solved to a specified order finite difference, 
//  i.e. either second or fourth order.
//
//  s(r, z) is a linear source.
//  f(r, z) is the RHS.
//
//  It returns the solution u(r, z) and the residual res(r, z).
//
extern "C" void flat_laplacian_(double *u,	// Output solution.
	const double *f,	// Input RHS.
	double *res,		// Ouput residual.
	const double *s,	// Input linear source.
	const double uInf,	// u value at infinity for Robin BC.
	const int robin,	// Robin BC type: 1, 2, 3.
	const int r_sym,	// R symmetry: 1(even), -1(odd).
	const int z_sym,	// Z symmetry: 1(even), -1(odd).
	const int NrInterior,	// Number of r interior points.
	const int NzInterior,	// Number of z interior points.
	const int ghost_zones,	// Number of ghost zones.
	const double dr, 	// Spatial step in r.
	const double dz,	// Spatial step in z.
	const int norder,	// Finite difference evolution: 2 or 4.
	const int precond_use) 	// Calculate and/or use preconditioner.
{
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
	double *g_s = (double *)malloc(g_size);
	double *g_res = (double *)malloc(g_size);

	// Reduce arrays.
	ghost_reduce(u, g_u, NrInterior, NzInterior, ghost);
	ghost_reduce(f, g_f, NrInterior, NzInterior, ghost);
	ghost_reduce(s, g_s, NrInterior, NzInterior, ghost);
	ghost_reduce(res, g_res, NrInterior, NzInterior, ghost);

	// Set new ghost.
	ghost = 1;

	// Set temporary total number of points.
	int NrTotal = NrInterior + 2;
	int NzTotal = NzInterior + 2;

	// Allocate and generate CSR matrix.
	csr_matrix A;
	int DIM0 = NrTotal * NzTotal;
	int nnz0 = nnz_flat_laplacian(NrInterior, NzInterior, norder, robin);
	csr_allocate(&A, DIM0, DIM0, nnz0);
	printf("FLAT LAPLACIAN: Allocated CSR matrix with %d rows, %d columns and %d nnz.\n", A.nrows, A.ncols, A.nnz);

	// Fill CSR matrix.
	csr_gen_flat_laplacian(A, NrInterior, NzInterior, norder, dr, dz, g_s, g_f, uInf, robin, r_sym, z_sym);

#ifdef DEBUG
	write_single_file(g_f, "f.asc", NrTotal, NzTotal);
#endif

	// Elliptic solver return variables.
	double norm = 0.0;
	int convergence = 0;
	double tol = (norder == 4) ? dr * dr * dz * dz : dr * dz;

	// Call elliptic solver.
	pardiso_wrapper(A, g_u, g_f, g_res, tol, &norm, &convergence, INFNORM, precond_use);

	// Check solver convergence.
	if (convergence == 1)
	{
		printf("FLAT LAPLACIAN: Solver converged!\n");
	}
	else
	{
		printf("FLAT LAPLACIAN: WARNING possible no convergence: %d.!\n", convergence);
	}
	printf("FLAT LAPLACIAN: ||r|| = %3.3E.\n", norm);

	// Reset ghost and total number of points.
	ghost = temp_ghost;
	NzTotal = NzInterior + ghost + 1;
	NrTotal = NrInterior + ghost + 1;

	// Transfer solution and residual to original arrays.
	ghost_fill(g_u, u, r_sym, z_sym, NrInterior, NzInterior, ghost);
	ghost_fill(g_res, res, r_sym, z_sym, NrInterior, NzInterior, ghost);

	// Clear memory.
	free(g_u);
	free(g_res);
	free(g_f);
	free(g_s);

	// Clear CSR matrix.
	csr_deallocate(&A);

	return;
}
