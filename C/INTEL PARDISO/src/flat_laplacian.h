void flat_laplacian(double *u,	// Output solution.
	double *res,		// Ouput residual.
	const double *s,	// Input linear source.
	const double *f,	// Input RHS.
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
	const int lr_use,	// Low rank update.
	const int perm_use,	// Calculate and/or use permutation.
	const int precond_use);	// Calculate and/or use preconditioner.
