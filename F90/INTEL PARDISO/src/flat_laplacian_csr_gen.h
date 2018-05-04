// Nonzero elements calculator.
int nnz_flat_laplacian(const int NrInterior, const int NzInterior, const int order, const int robin);

// Write CSR matrix for the flat laplacian.
void csr_gen_flat_laplacian(csr_matrix A,
	const int NrInterior,		// Number of r interior points.
	const int NzInterior,		// Number of z interior points.
	const int order,		// Finite difference order: 2 or 4.
	const double dr,		// Spatial step in r.
	const double dz,		// Spatial step in z.
	const double *s,		// Linear source.
	double *f,			// RHS.
	const double uInf,		// Value at infinity.
	const int robin,		// Robin BC type: 1, 2, 3.
	const int r_sym,		// R symmetry: 1(even), -1(odd).
	const int z_sym);		// Z symmetry: 1(even), -1(odd).
