// General elliptic equation, solves the linear equation:
//     2       2       2          
// (a d  +  b d  +  c d  +  d d  +  e d  +  s) u = f,
//     rr      rz      zz      r       z
//
// where u is the solution and a, b, c, d, e, f, s are all 
// functions of (r, z).
// 
void general_elliptic(double *u,// Output solution.
	double *res,                // Output residual. 
	const double *ell_a,        // Input a coefficient.
	const double *ell_b,        // Input b coefficient.
	const double *ell_c,        // Input c coefficient.
	const double *ell_d,        // Input d coefficient.
	const double *ell_e,        // Input e coefficient.
	const double *ell_s,        // Input s coefficient.
	const double *ell_f,        // Input f coefficient.
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
	const int lr_use,	// Use low rank update.
	const int precond_use); 	// Calculate and/or use preconditioner.
