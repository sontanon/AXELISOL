// Non-zero calculator.
int nnz_general_elliptic(const int NrInterior, const int NzInterior, const int order, const int robin);

// Write CSR matrix for a general elliptic equation.
void csr_gen_general_elliptic(csr_matrix A,	// CSR matrix structure.
	const int NrInterior,			// Number of r interior points.	
	const int NzInterior,			// Number of z interior points.
	const int order,			    // Finite difference order: 2 or 4.
	const double dr,			    // Spatial step in r.
	const double dz,			    // Spatial step in z.
	const double *ell_a,			// Coefficient of (d^2/dr^2)
	const double *ell_b,			// Coefficient of (d^2/drdz)
	const double *ell_c,			// Coefficient of (d^2/dz^)
	const double *ell_d,			// Coefficient of (d/dr)
	const double *ell_e,			// Coefficient of (d/dz)
	const double *ell_s,			// Linear source.
	double *ell_f,				    // RHS.
	const double uInf,			    // Value at infinity.
	const int robin,			    // Robin BC type: 1, 2, 3.
	const int r_sym,			    // R symmetry: 1(even), -1(odd).
	const int z_sym);			    // Z symmetry: 1(even), -1(odd).
