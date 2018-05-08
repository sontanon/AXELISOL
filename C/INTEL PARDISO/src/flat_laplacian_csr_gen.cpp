// Global header.
// One-based indexing BASE is defined in this header.
#include "tools.h"

// Print CSR matrix for debug.
#undef DEBUG

// Nonzero elements calculator.
int nnz_flat_laplacian(const int NrInterior, const int NzInterior, const int order, const int robin)
{
	int nnz;

	int n_robin = robin + order;

	// Second order Laplacian.
	if (order == 2)
	{
		nnz = 5 * NrInterior * NzInterior
			+ (2 + n_robin) * (NrInterior + NzInterior)
			+ (2 + 2 + 2 + n_robin);
	}
	// Fourth-order Laplacian.
	else if (order == 4)
	{
		nnz = 9 * (NrInterior - 2) * (NzInterior - 2)
			+ (8 + 10) * (NrInterior + NzInterior - 4)
			+ (7 + 9 + 9 + 11)
			+ (2 + n_robin) * (NrInterior + NzInterior)
			+ (2 + 2 + 2 + n_robin);
	}
	return nnz;
}

// Write CSR matrix for the flat laplacian.
void csr_gen_flat_laplacian(csr_matrix A, // CSR matrix structure.
	const int NrInterior,		// Number of r interior points.
	const int NzInterior,		// Number of z interior points.
	const int order,		    // Finite difference order: 2 or 4.
	const double dr,		    // Spatial step in r.
	const double dz,		    // Spatial step in z.
	const double *s,		    // Linear source.
	double *f,			        // RHS.
	const double uInf,          // Value at infinity.
	const int robin,		    // Robin BC type: 1, 2, 3.
	const int r_sym,		    // R symmetry: 1(even), -1(odd).
	const int z_sym)		    // Z symmetry: 1(even), -1(odd).
{
	// Constant numbers.
	const double third = 1.0 / 3.0;
	const double sixth = 1.0 / 6.0;
	const double twelfth = 1.0 / 12.0;

	// Grid extensions.
	int NrTotal = NrInterior + 2;
	int NzTotal = NzInterior + 2;

	// Number of nonzero elements.
	int nnz = nnz_flat_laplacian(NrInterior, NzInterior, order, robin);

	// Number of Robin elements.
	int n_robin = robin + order;

	// Number of elements we have filled in.
	int offset = 0;

	// Auxiliary variables.
	double r, z, rr2, ir, rrodrr;
	double robin1, robin2, robin3;
	int i, j, t_offset;

	// Ratios for variable step size.
	double roz = dr / dz;
	double zor = dz / dr;

	// SECOND-ORDER FLAT LAPLACIAN.
	if (order == 2)
	{
		// Lower-left corner: diagonal symmetry.
		A.ia[IDX(0, 0)] = BASE + offset;
		A.a[offset] = 1.0;
		A.a[offset + 1] = -(double)(r_sym * z_sym);
		A.ja[offset] = BASE + IDX(0, 0);
		A.ja[offset + 1] = BASE + IDX(1, 1);
		f[IDX(0, 0)] = 0.0;
		offset += 2;

		// Set temporary offset.
		t_offset = offset;

		// Fill left-boundary using axis symmetry.
		#pragma omp parallel shared(A, f) private(offset)
		{
			#pragma omp for schedule(guided)
			for (j = 1; j < NzInterior + 1; j++)
			{
				// Each j iteration fills 2 elements.
				offset = t_offset + 2 * (j - 1);
				A.ia[IDX(0, j)] = BASE + offset;
				A.a[offset] = 1.0;
				A.a[offset + 1] = -(double)r_sym;
				A.ja[offset] = BASE + IDX(0, j);
				A.ja[offset + 1] = BASE + IDX(1, j);
				f[IDX(0, j)] = 0.0;
			}
		}

		// We have now filled:
		offset = 2 + 2 * NzInterior;

		// Upper-left corner: also axis symmetry.
		A.ia[IDX(0, NzInterior + 1)] = BASE + offset;
		A.a[offset] = 1.0;
		A.a[offset + 1] = -(double)r_sym;
		A.ja[offset] = BASE + IDX(0, NzInterior + 1);
		A.ja[offset + 1] = BASE + IDX(1, NzInterior + 1);
		f[IDX(0, NzInterior + 1)] = 0.0;
		offset += 2;

		// Set temporary offset.
		t_offset = offset;

		// Now come the interior points plus the top and bottom boundaries with
		// Robin and equatorial symmetry respectively.
		#pragma omp parallel shared(A, f) private(offset, j, r, z, ir, rr2, robin1, robin2, robin3)
		{
			#pragma omp for schedule(guided)
			for (i = 1; i < NrInterior + 1; i++)
			{
				// Each iteration of i loop will fill 5 * NzInterior + (2 + n_robin) values.
				offset = t_offset + (i - 1) * (5 * NzInterior + 2 + n_robin);

				// R coordinate.
				r = (double)i - 0.5;
				// Inverse.
				ir = 1.0 / r;

				// Do bottom boundary first with equatorial symmetry.
				A.ia[IDX(i, 0)] = BASE + offset;
				A.a[offset] = 1.0;
				A.a[offset + 1] = -(double)z_sym;
				A.ja[offset] = BASE + IDX(i, 0);
				A.ja[offset + 1] = BASE + IDX(i, 1);
				f[IDX(i, 0)] = 0.0;
				offset += 2;

				// Now loop over interior points.
				for (j = 1; j < NzInterior + 1; j++)
				{
					// Row begins at offset.
					A.ia[IDX(i, j)] = BASE + offset;
					// Values.
					A.a[offset] = (1.0 - 0.5 * ir) * zor;
					A.a[offset + 1] = 1.0 * roz;
					A.a[offset + 2] = dr * dz * s[IDX(i, j)] - 2.0 * (roz + zor);
					A.a[offset + 3] = 1.0 * roz;
					A.a[offset + 4] = (1.0 + 0.5 * ir) * zor;
					// Columns.
					A.ja[offset] = BASE + IDX(i - 1, j);
					A.ja[offset + 1] = BASE + IDX(i, j - 1);
					A.ja[offset + 2] = BASE + IDX(i, j);
					A.ja[offset + 3] = BASE + IDX(i, j + 1);
					A.ja[offset + 4] = BASE + IDX(i + 1, j);
					// Multiply RHS by dr * dz.
					f[IDX(i, j)] *= dr * dz;
					// Increase offset by five.
					offset += 5;
				}

				// Now fill the top boundary with Robin.
				j = NzInterior + 1;
				// Z coordinate.
				z = (double)j - 0.5;
				// Radial coordinate.
				rr2 = r * r * roz * roz + z * z;
				A.ia[IDX(i, NzInterior + 1)] = BASE + offset;
				switch (robin)
				{
					case 1:
						A.a[offset] = 0.5 * rr2 / z;
						A.a[offset + 1] = -2.0 * rr2 / z;
						A.a[offset + 2] = 1.0 + 1.5 * rr2 / z;
						A.ja[offset] = BASE + IDX(i, NzInterior - 1);
						A.ja[offset + 1] = BASE + IDX(i, NzInterior);
						A.ja[offset + 2] = BASE + IDX(i, NzInterior + 1);
						break;
					case 2:
						robin2 = (rr2 / z) * (rr2 / z);
						robin1 = (rr2 / z) * (4.0 - (r * roz / z) * (r * roz / z));
						A.a[offset] = -0.5 * robin2;
						A.a[offset + 1] = 2.0 * robin2 + 0.25 * robin1;
						A.a[offset + 2] = -2.5 * robin2 - robin1;
						A.a[offset + 3] = 1.0 + robin2 + 0.75 * robin1;
						A.ja[offset] = BASE + IDX(i, NzInterior - 2);
						A.ja[offset + 1] = BASE + IDX(i, NzInterior - 1);
						A.ja[offset + 2] = BASE + IDX(i, NzInterior);
						A.ja[offset + 3] = BASE + IDX(i, NzInterior + 1);
						break;
					case 3:
						robin3 = (rr2 / z) * (rr2 / z) * (rr2 / z);
						robin2 = (rr2 / z) * (rr2 / z) * (9.0 - 3.0 * (r * roz / z) * (r * roz / z));
						robin1 = (rr2 / z) * (18.0 + (r * roz / z) * (r * roz / z) * (-9.0 + 3.0 * (rr2 / (z * z))));
						A.a[offset] = 0.25 * robin3;
						A.a[offset + 1] = -(7.0 * robin3 + robin2) / 6.0;
						A.a[offset + 2] = 2.0 * robin3 + 2.0 * robin2 / 3.0 + robin1 / 12.0;
						A.a[offset + 3] = -(1.5 * robin3 + 5.0 * robin2 / 6.0 + robin1 / 3.0);
						A.a[offset + 4] = 1.0 + 5.0 * robin3 / 12.0 + robin2 / 3.0 + 0.25 * robin1;
						A.ja[offset] = BASE + IDX(i, NzInterior - 3);
						A.ja[offset + 1] = BASE + IDX(i, NzInterior - 2);
						A.ja[offset + 2] = BASE + IDX(i, NzInterior - 1);
						A.ja[offset + 3] = BASE + IDX(i, NzInterior);
						A.ja[offset + 4] = BASE + IDX(i, NzInterior + 1);
						break;
				}
				// Also fill RHS term.
				f[IDX(i, NzInterior + 1)] = uInf;
			}
		}

		// At this point we have now filled:
		offset = 4 + 2 * NzInterior + 5 * NrInterior * NzInterior + (2 + n_robin) * NrInterior;

		// Lower-right corner: equatorial symmetry.
		A.ia[IDX(NrInterior + 1, 0)] = BASE + offset;
		A.a[offset] = 1.0;
		A.a[offset + 1] = -(double)z_sym;
		A.ja[offset] = BASE + IDX(NrInterior + 1, 0);
		A.ja[offset + 1] = BASE + IDX(NrInterior + 1, 1);
		f[IDX(NrInterior + 1, 0)] = 0.0;
		offset += 2;

		// Set temporary offset.
		t_offset = offset;

		// Robin boundary.
		r = (double)NrInterior + 0.5;
		#pragma omp parallel shared(A, f) private(offset, z, rr2, robin1, robin2, robin3)
		{
			#pragma omp for schedule(guided)
			for (j = 1; j < NzInterior + 1; j++)
			{

				z = (double)j - 0.5;
				rr2 = r * r + z * z * zor * zor;
				// Each iteration of the loop fills n_robin elements.
				offset = t_offset + n_robin * (j - 1);

				A.ia[IDX(NrInterior + 1, j)] = BASE + offset;
				switch (robin)
				{
					case 1:
						A.a[offset] = 0.5 * rr2 / r;
						A.a[offset + 1] = -2.0 * rr2 / r;
						A.a[offset + 2] = 1.0 + 1.5 * rr2 / r;
						A.ja[offset] = BASE + IDX(NrInterior - 1, j);
						A.ja[offset + 1] = BASE + IDX(NrInterior, j);
						A.ja[offset + 2] = BASE + IDX(NrInterior + 1, j);
						break;
					case 2:
						robin2 = (rr2 / r) * (rr2 / r);
						robin1 = (rr2 / r) * (4.0 - (z * zor / r) * (z * zor / r));
						A.a[offset] = -0.5 * robin2;
						A.a[offset + 1] = 2.0 * robin2 + 0.25 * robin1;
						A.a[offset + 2] = -2.5 * robin2 - robin1;
						A.a[offset + 3] = 1.0 + robin2 + 0.75 * robin1;
						A.ja[offset] = BASE + IDX(NrInterior - 2, j);
						A.ja[offset + 1] = BASE + IDX(NrInterior - 1, j);
						A.ja[offset + 2] = BASE + IDX(NrInterior, j);
						A.ja[offset + 3] = BASE + IDX(NrInterior + 1, j);
						break;
					case 3:
						robin3 = (rr2 / r) * (rr2 / r) * (rr2 / r);
						robin2 = (rr2 / r) * (rr2 / r) * (9.0 - 3.0 * (z * zor / r) * (z * zor / r));
						robin1 = (rr2 / r) * (18.0 + (z * zor / r) * (z * zor / r) * (-9.0 + 3.0 * (rr2 / (r * r))));
						A.a[offset] = 0.25 * robin3;
						A.a[offset + 1] = -(7.0 * robin3 + robin2) / 6.0;
						A.a[offset + 2] = 2.0 * robin3 + 2.0 * robin2 / 3.0 + robin1 / 12.0;
						A.a[offset + 3] = -(1.5 * robin3 + 5.0 * robin2 / 6.0 + robin1 / 3.0);
						A.a[offset + 4] = 1.0 + 5.0 * robin3 / 12.0 + robin2 / 3.0 + 0.25 * robin1;
						A.ja[offset] = BASE + IDX(NrInterior - 3, j);
						A.ja[offset + 1] = BASE + IDX(NrInterior - 2, j);
						A.ja[offset + 2] = BASE + IDX(NrInterior - 1, j);
						A.ja[offset + 3] = BASE + IDX(NrInterior, j);
						A.ja[offset + 4] = BASE + IDX(NrInterior + 1, j);
						break;
				}
				// Also fill RHS term.
				f[IDX(NrInterior + 1, j)] = uInf;
			}
		}

		// At this point, we have now filled:
		offset = 6 + (2 + n_robin) * (NzInterior + NrInterior) + 5 * NrInterior * NzInterior;

		// Upper-right corner: fill with Robin.
		r = (double)NrInterior + 0.5;
		z = (double)NzInterior + 0.5;
		rrodrr = sqrt((r * r * dr * dr + z * z * dz * dz) / (dr * dr + dz * dz));
		A.ia[IDX(NrInterior + 1, NzInterior + 1)] = BASE + offset;
		switch (robin)
		{
			case 1:
				A.a[offset] = 0.5 * rrodrr;
				A.a[offset + 1] = -2.0 * rrodrr;
				A.a[offset + 2] = 1.0 + 1.5 * rrodrr;
				A.ja[offset] = BASE + IDX(NrInterior - 1, NzInterior - 1);
				A.ja[offset + 1] = BASE + IDX(NrInterior, NzInterior);
				A.ja[offset + 2] = BASE + IDX(NrInterior + 1, NzInterior + 1);
				break;
			case 2:
				robin2 = rrodrr * rrodrr;
				robin1 = 4.0 * rrodrr;
				A.a[offset] = -0.5 * robin2;
				A.a[offset + 1] = 2.0 * robin2 + 0.25 * robin1;
				A.a[offset + 2] = -2.5 * robin2 - robin1;
				A.a[offset + 3] = 1.0 + robin2 + 0.75 * robin1;
				A.ja[offset] = BASE + IDX(NrInterior - 2, NzInterior - 2);
				A.ja[offset + 1] = BASE + IDX(NrInterior - 1, NzInterior -1);
				A.ja[offset + 2] = BASE + IDX(NrInterior, NzInterior);
				A.ja[offset + 3] = BASE + IDX(NrInterior + 1, NzInterior + 1);
				break;
			case 3:
				robin3 = rrodrr * rrodrr * rrodrr;
				robin2 = 9.0 * rrodrr * rrodrr;
				robin1 = 18.0 * rrodrr;
				A.a[offset] = 0.25 * robin3;
				A.a[offset + 1] = -(7.0 * robin3 + robin2) / 6.0;
				A.a[offset + 2] = 2.0 * robin3 + 2.0 * robin2 / 3.0 + robin1 / 12.0;
				A.a[offset + 3] = -(1.5 * robin3 + 5.0 * robin2 / 6.0 + robin1 / 3.0);
				A.a[offset + 4] = 1.0 + 5.0 * robin3 / 12.0 + robin2 / 3.0 + 0.25 * robin1;
				A.ja[offset] = BASE + IDX(NrInterior - 3, NzInterior - 3);
				A.ja[offset + 1] = BASE + IDX(NrInterior - 2, NzInterior - 2);
				A.ja[offset + 2] = BASE + IDX(NrInterior - 1, NzInterior -1);
				A.ja[offset + 3] = BASE + IDX(NrInterior, NzInterior);
				A.ja[offset + 4] = BASE + IDX(NrInterior + 1, NzInterior + 1);
				break;
		}
		offset += n_robin;
		// Also fill RHS term.
		f[IDX(NrInterior + 1, NzInterior + 1)] = uInf;

		// Assert fill-in.
		assert(nnz == offset);

		// Fill last element of row offsets.
		A.ia[IDX(NrInterior + 1, NzInterior + 1) + 1] = BASE + nnz;
	}
	// FOURTH-ORDER FLAT LAPLACIAN.
	else if (order == 4)
	{
		// Set offset to zero.
		offset = 0;

		// Lower-left corner: diagonal symmetry.
		A.ia[IDX(0, 0)] = BASE + offset;
		A.a[offset] = 1.0;
		A.a[offset + 1] = -(double)(r_sym * z_sym);
		A.ja[offset] = BASE + IDX(0, 0);
		A.ja[offset + 1] = BASE + IDX(1, 1);
		f[IDX(0, 0)] = 0.0;
		offset += 2;

		// offset = 2.

		// Set temporary offset.
		t_offset = offset;

		// Fill left-boundary using axis symmetry.
		#pragma omp parallel shared(A, f) private(offset)
		{
			#pragma omp for schedule(guided)
			for (j = 1; j < NzInterior + 1; j++)
			{
				// Each j iteration fills 2 elements.
				offset = t_offset + 2 * (j - 1);
				A.ia[IDX(0, j)] = BASE + offset;
				A.a[offset] = 1.0;
				A.a[offset + 1] = -(double)r_sym;
				A.ja[offset] = BASE + IDX(0, j);
				A.ja[offset + 1] = BASE + IDX(1, j);
				f[IDX(0, j)] = 0.0;
			}
		}

		// We have now filled:
		offset = 2 + 2 * NzInterior;

		// offset = 2 + NzInterior * 2.

		// Upper-left corner: also axis symmetry.
		A.ia[IDX(0, NzInterior + 1)] = BASE + offset;
		A.a[offset] = 1.0;
		A.a[offset + 1] = -(double)r_sym;
		A.ja[offset] = BASE + IDX(0, NzInterior + 1);
		A.ja[offset + 1] = BASE + IDX(1, NzInterior + 1);
		f[IDX(0, NzInterior + 1)] = 0.0;
		offset += 2;

		// offset = (2 + 2) + 2 * NzInterior.

		// Left-interior points with semi-one sided finite difference.
		i = 1;

		// Bottom boundary with equatorial symmetry.
		j = 0;
		A.ia[IDX(1, 0)] = BASE + offset;
		A.a[offset] = 1.0;
		A.a[offset + 1] = -(double)z_sym;
		A.ja[offset] = BASE + IDX(1, 0);
		A.ja[offset + 1] = BASE + IDX(1, 1);
		f[IDX(1, 0)] = 0.0;
		offset += 2;

		// offset = (2 + 2) + 2 * NzInterior 
		// + 2.

		// Semi-one sided bottom with extra ghost zone.
		//
		//    o
		//    |
		//    o
		//    |
		// o--x--o--o
		//    |
		//    o
		//
		// Thus:
		// 
		// u(i-1, j  ) : 16/12 - 8h/12r.
		// u(i  , j-1) : 16/12
		// u(i  , j  ) : s(i, j) * h**2 - 30/12 - 30/12
		// u(i  , j+1) : 16/12 + (z_sym)*(-1/12)
		// u(i  , j+2) : -1/12
		// u(i+1, j  ) : 16/12 + 8h/12r + (r_sym)*(-1/12 + 1h/12r)
		// u(i+2, j  ) : -1/12 - 1h/12r
		//
		j = 1;
		r = 0.5;
		ir = 2.0;

		A.ia[IDX(1, 1)] = BASE + offset;
		A.a[offset] = (4.0 * third - 2.0 * third * ir) * zor;
		A.a[offset + 1] = 4.0 * third * roz;
		A.a[offset + 2] = dr * dz * s[IDX(1, 1)] - 2.5 * (zor + roz);
		A.a[offset + 3] = (4.0 * third + (double)z_sym * (-twelfth)) * roz;
		A.a[offset + 4] = -twelfth * roz;
		A.a[offset + 5] = (4.0 * third + 2.0 * third * ir + (double)r_sym * (-twelfth + twelfth * ir)) * zor;
		A.a[offset + 6] = (-twelfth - twelfth * ir) * zor;
		A.ja[offset] = BASE + IDX(0, 1);
		A.ja[offset + 1] = BASE + IDX(1, 0);
		A.ja[offset + 2] = BASE + IDX(1, 1);
		A.ja[offset + 3] = BASE + IDX(1, 2);
		A.ja[offset + 4] = BASE + IDX(1, 3);
		A.ja[offset + 5] = BASE + IDX(2, 1);
		A.ja[offset + 6] = BASE + IDX(3, 1);
		f[IDX(1, 1)] *= dr * dz;
		offset += 7;

		// offset = (2 + 2) +  2 * NzInterior 
		// + (2 + 7).

		// Set temporary offset.
		t_offset = offset;

		// Combination of semi-onesided and centered.
		//
		//    o
		//    |
		//    o
		//    |
		// o--x--o--o
		//    |
		//    o
		//    |
		//    o
		//
		// Thus:
		//
		// u(i-1, j  ) : 16/12 - 8h/12r
		// u(i  , j-2) : -1/12
		// u(i  , j-1) : 16/12
		// u(i  , j  ) : s(i, j) * h**2 - 30/12 - 30/12
		// u(i  , j+1) : 16/12
		// u(i  , j+2) : -1/12
		// u(i+1, j  ) : 16/12 + 8h/12r + (r_sym)*(-1/12 + 1h/12r)
		// u(i+2, j  ) : -1/12 - 1h/12r
		//
		#pragma omp parallel shared(A, f) private(offset)
		{
			#pragma omp for schedule(guided)
			for (j = 2; j < NzInterior; j++)
			{
				// Each iteration fills 8 elements.
				offset = t_offset + 8 * (j - 2);
				A.ia[IDX(1, j)] = BASE + offset;
				A.a[offset] = (4.0 * third - 2.0 * third * ir) * zor;
				A.a[offset + 1] = -twelfth * roz;
				A.a[offset + 2] = 4.0 * third * roz;
				A.a[offset + 3] = dr * dz * s[IDX(1, j)] - 2.5 * (zor + roz);
				A.a[offset + 4] = 4.0 * third * roz;
				A.a[offset + 5] = -twelfth * roz; 
				A.a[offset + 6] = (4.0 * third + 2.0 * third * ir + (double)r_sym * (-twelfth + twelfth * ir)) * zor;
				A.a[offset + 7] = (-twelfth - twelfth * ir) * zor;
				A.ja[offset] = BASE + IDX(0, j);
				A.ja[offset + 1] = BASE + IDX(1, j - 2);
				A.ja[offset + 2] = BASE + IDX(1, j - 1);
				A.ja[offset + 3] = BASE + IDX(1, j);
				A.ja[offset + 4] = BASE + IDX(1, j + 1);
				A.ja[offset + 5] = BASE + IDX(1, j + 2);
				A.ja[offset + 6] = BASE + IDX(2, j);
				A.ja[offset + 7] = BASE + IDX(3, j);
				f[IDX(1, j)] *= dr * dz;
			}
		}

		// At this point we have filled:
		offset = (2 + 2) + 2 * NzInterior
			+ 2 + 7 + 8 * (NzInterior - 2);

		// offset = (2 + 2) + 2 * NzInterior 
		// + 2 + 7 + 8 * (NzInterior - 2).


		// Semi-one sided top.
		//
		//    o
		//    |
		// o--x--o--o
		//    |
		//    o
		//    |
		//    o
		//    |
		//    o
		//    |
		//    o
		//
		// Recall:
		// - Semi-one sided on z.
		// u''(x) = (10u(i + 1) - 15u(i) - 4u(i - 1) + 14u(i - 2) - 6u(i - 3) + u(i - 4)/12h**2
		// 
		// Thus:
		//
		// u(i-1, j  ) : 16/12 - 8h/12r
		// u(i  , j-4) : 1/12
		// u(i  , j-3) : -6/12
		// u(i  , j-2) : 14/12
		// u(i  , j-1) : -4/12
		// u(i  , j  ) : s(i, j) * h **2  - 30/12 - 15/12
		// u(i  , j+1) : 10/12
		// u(i+1, j  ) : 16/12 + 8h/12r + (r_sym)*(-1/12 + 1h/12r)
		// u(i+2, j  ) : -1/12 - 1h/12r
		//
		j = NzInterior;
		A.ia[IDX(1, NzInterior)] = BASE + offset;
		A.a[offset] = (4.0 * third - 2.0 * third * ir) * zor;
		A.a[offset + 1] = twelfth * roz;
		A.a[offset + 2] = -0.5 * roz;
		A.a[offset + 3] = 7.0 * sixth * roz;
		A.a[offset + 4] = -third * roz;
		A.a[offset + 5] = dr * dz * s[IDX(i, j)] - 2.5 * zor - 1.25 * roz;
		A.a[offset + 6] = 5.0 * sixth * roz;
		A.a[offset + 7] = (4.0 * third + 2.0 * third * ir + (double)r_sym * (-twelfth + twelfth * ir)) * zor;
		A.a[offset + 8] = (-twelfth - twelfth * ir) * zor;
		A.ja[offset] = BASE + IDX(0, NzInterior);
		A.ja[offset + 1] = BASE + IDX(1, NzInterior - 4);
		A.ja[offset + 2] = BASE + IDX(1, NzInterior - 3);
		A.ja[offset + 3] = BASE + IDX(1, NzInterior - 2);
		A.ja[offset + 4] = BASE + IDX(1, NzInterior - 1);
		A.ja[offset + 5] = BASE + IDX(1, NzInterior);
		A.ja[offset + 6] = BASE + IDX(1, NzInterior + 1);
		A.ja[offset + 7] = BASE + IDX(2, NzInterior);
		A.ja[offset + 8] = BASE + IDX(3, NzInterior);
		f[IDX(1, NzInterior)] *= dr * dz;
		offset += 9;

		// offset = (2 + 2) + 2 * NzInterior 
		// + 2 + 7 + 8 * (NzInterior - 2) + 9.

		j = NzInterior + 1;
		z = (double)j - 0.5;
		rr2 = r * r * roz * roz + z * z;
		A.ia[IDX(1, NzInterior + 1)] = BASE + offset;
		switch (robin)
		{
			case 1:
				robin1 = rr2 / z;
				A.a[offset] = robin1 / 4.0;
				A.a[offset + 1] = -4.0 * robin1 / 3.0;
				A.a[offset + 2] = 3.0 * robin1;
				A.a[offset + 3] = -4.0 * robin1;
				A.a[offset + 4] = 1.0 + 25.0 * robin1 / 12.0;
				A.ja[offset] = BASE + IDX(i, NzInterior - 3);
				A.ja[offset + 1] = BASE + IDX(i, NzInterior - 2);
				A.ja[offset + 2] = BASE + IDX(i, NzInterior - 1);
				A.ja[offset + 3] = BASE + IDX(i, NzInterior);
				A.ja[offset + 4] = BASE + IDX(i, NzInterior + 1);
				break;
		}
		f[IDX(1, NzInterior + 1)] = uInf;
		offset += n_robin;

		// offset = (2 + 2) + 2 * NzInterior 
		// + 2 + 7 + 8 * (NzInterior - 2) + 9 + n_robin.

		// Set temporary offset.
		t_offset = offset;

		#pragma omp parallel shared(A, f) private(j, offset, r, z, ir, rr2,\
		robin1, robin2, robin3)
		{
			#pragma omp for schedule(guided)
			for (i = 2; i < NrInterior; i++)
			{
				// Each iteration fills 2 + 8 + 9 * (NzInterior - 2) + 10 + n_robin points.
				offset = t_offset + (i - 2) * (9 * NzInterior + 2 + n_robin);

				// R coordinate.
				r = (double)i - 0.5;

				// Inverse.
				ir = 1.0 / r;

				// Equatorial symmetry.
				j = 0;
				A.ia[IDX(i, 0)] = BASE + offset;
				A.a[offset] = 1.0;
				A.a[offset + 1] = -(double)z_sym;
				A.ja[offset] = BASE + IDX(i, 0);
				A.ja[offset + 1] = BASE + IDX(i, 1);
				f[IDX(i, 0)] = 0.0;
				offset += 2;

				// Semi-onesided and centered difference.
				//
				//       o
				//       |
				//       o
				//       |
				// o--o--x--o--o
				//       |
				//       o
				//
				// Thus:
				// 
				// u(i-2, j  ) : -1/12 + 1h/12r
				// u(i-1, j  ) : 16/12 - 8h/12r
				// u(i  , j-1) : 16/12
				// u(i  , j  ) : s(i,j) * h**2 -30/12 - 30/12
				// u(i  , j+1) : 16/12 + (z_sym)*(-1/12)
				// u(i  , j+2) : -1/12
				// u(i+1, j  ) : 16/12 + 8h/12r
				// u(i+2, j  ) : -1/12 - 1h/12r
				//
				j = 1;
				A.ia[IDX(i, j)] = BASE + offset;
				A.a[offset] = (-twelfth + twelfth * ir) * zor;
				A.a[offset + 1] = (4.0 * third - 2.0 * third * ir) * zor;
				A.a[offset + 2] = 4.0 * third * roz;
				A.a[offset + 3] = dr * dz * s[IDX(i, j)] - 2.5 * (roz + zor);
				A.a[offset + 4] = (4.0 * third + (double)z_sym * (-twelfth)) * roz;
				A.a[offset + 5] = -twelfth * roz;
				A.a[offset + 6] = (4.0 * third + 2.0 * third * ir) * zor;
				A.a[offset + 7] = (-twelfth - twelfth * ir) * zor;
				A.ja[offset] = BASE + IDX(i - 2, 1);
				A.ja[offset + 1] = BASE + IDX(i - 1, 1);
				A.ja[offset + 2] = BASE + IDX(i, 0);
				A.ja[offset + 3] = BASE + IDX(i, 1);
				A.ja[offset + 4] = BASE + IDX(i, 2);
				A.ja[offset + 5] = BASE + IDX(i, 3);
				A.ja[offset + 6] = BASE + IDX(i + 1, 1);
				A.ja[offset + 7] = BASE + IDX(i + 2, 1);
				f[IDX(i, j)] *= dr * dz;
				offset += 8;

				// Centered interior points.
				//
				//       o
				//       |
				//       o
				//       |
				// o--o--x--o--o
				//       |
				//       o
				//       |
				//       o
				//
				// Recall:
				// - Centered on r:
				// u' (x) = (-u(i+2)+8u(i+1)-8u(i-1)+u(i-2))/12h
				// u''(x) = (-(u(i+2)+u(i-2))+16(u(i+1)+u(i-1))-30u(i))/12h**2
				// - Centered on z:
				// u''(x) = (-(u(i+2)+u(i-2))+16(u(i+1)+u(i-1))-30u(i))/12h**2
				//
				// Thus:
				//
				// u(i-2, j  ) : -1/12 + 1h/12r
				// u(i-1, j  ) : 16/12 - 8h/12r
				// u(i  , j-2) : -1/12
				// u(i  , j-1) : 16/12
				// u(i  , j  ) : s(i,j) * h**2 -30/12 - 30/12
				// u(i  , j+1) : 16/12
				// u(i  , j+2) : -1/12
				// u(i+1, j  ) : 16/12 + 8h/12r
				// u(i+2, j  ) : -1/12 - 1h/12r
				//
				for (j = 2; j < NzInterior; j++)
				{
					A.ia[IDX(i, j)] = BASE + offset;
					A.a[offset] = (-twelfth + twelfth * ir) * zor;
					A.a[offset + 1] = (4.0 * third - 2.0 * third * ir) * zor;
					A.a[offset + 2] = -twelfth * roz;
					A.a[offset + 3] = 4.0 * third * roz;
					A.a[offset + 4] = dr * dz * s[IDX(i, j)] - 2.5 * (zor + roz);
					A.a[offset + 5] = 4.0 * third * roz;
					A.a[offset + 6] = -twelfth * roz;
					A.a[offset + 7] = (4.0 * third + 2.0 * third * ir) * zor;
					A.a[offset + 8] = (-twelfth - twelfth * ir) * zor;
					A.ja[offset] = BASE + IDX(i - 2, j);
					A.ja[offset + 1] = BASE + IDX(i - 1, j);
					A.ja[offset + 2] = BASE + IDX(i, j - 2);
					A.ja[offset + 3] = BASE + IDX(i, j - 1);
					A.ja[offset + 4] = BASE + IDX(i, j);
					A.ja[offset + 5] = BASE + IDX(i, j + 1);
					A.ja[offset + 6] = BASE + IDX(i, j + 2);
					A.ja[offset + 7] = BASE + IDX(i + 1, j);
					A.ja[offset + 8] = BASE + IDX(i + 2, j);
					f[IDX(i, j)] *= dr * dz;
					offset += 9;
				}

				// Semi-onesided and centered difference.
				// 
				//       o
				//       |
				// o--o--x--o--o
				//       |
				//       o
				//       |
				//       o
				//       |
				//       o
				//       |
				//       o
				//
				// - Centered on r:
				// u' (x) = (-u(i+2)+8u(i+1)-8u(i-1)+u(i-2))/12h
				// u''(x) = (-(u(i+2)+u(i-2))+16(u(i+1)+u(i-1))-30u(i))/12h**2
				// - Semi-onesided on z:
				// u''(x) = (10u(i + 1) - 15u(i) - 4u(i - 1) + 14u(i - 2) - 6u(i - 3) + u(i - 4)/12h**2
				//
				// Thus:
				//
				// u(i-2, j  ) : -1/12 + 1h/12r
				// u(i-1, j  ) : 16/12 - 8h/12r
				// u(i  , j-4) : 1/12
				// u(i  , j-3) : -6/12
				// u(i  , j-2) : 14/12
				// u(i  , j-1) : -4/12
				// u(i  , j  ) : s(i,j) * h**2 -30/12 -15/12
				// u(i  , j+1) : 10/12
				// u(i+1, j  ) : 16/12 + 8h/12r
				// u(i+2, j  ) : -1/12 - 1h/12r
				//
				j = NzInterior;
				A.ia[IDX(i, j)] = BASE + offset;
				A.a[offset] = (-twelfth + twelfth * ir) * zor;
				A.a[offset + 1] = (4.0 * third - 2.0 * third * ir) * zor;
				A.a[offset + 2] = twelfth * roz;
				A.a[offset + 3] = -0.5 * roz;
				A.a[offset + 4] = 7.0 * sixth * roz;
				A.a[offset + 5] = -third * roz;
				A.a[offset + 6] = dr * dz * s[IDX(i, j)] - 2.5 * zor - 1.25 * roz;
				A.a[offset + 7] = 5.0 * sixth * roz;
				A.a[offset + 8] = (4.0 * third + 2.0 * third * ir) * zor;
				A.a[offset + 9] = (-twelfth - twelfth * ir) * zor;
				A.ja[offset] = BASE + IDX(i - 2, j);
				A.ja[offset + 1] = BASE + IDX(i - 1, j);
				A.ja[offset + 2] = BASE + IDX(i, j - 4);
				A.ja[offset + 3] = BASE + IDX(i, j - 3);
				A.ja[offset + 4] = BASE + IDX(i, j - 2);
				A.ja[offset + 5] = BASE + IDX(i, j - 1);
				A.ja[offset + 6] = BASE + IDX(i, j);
				A.ja[offset + 7] = BASE + IDX(i, j + 1);
				A.ja[offset + 8] = BASE + IDX(i + 1, j);
				A.ja[offset + 9] = BASE + IDX(i + 2, j);
				f[IDX(i, j)] *= dr * dz;
				offset += 10;

				j = NzInterior + 1;
				z = (double)j - 0.5;
				rr2 = r * r * roz * roz + z * z;
				A.ia[IDX(i, j)] = BASE + offset;
				switch (robin)
				{
				case 1:
					robin1 = rr2 / z;
					A.a[offset] = robin1 / 4.0;
					A.a[offset + 1] = -4.0 * robin1 / 3.0;
					A.a[offset + 2] = 3.0 * robin1;
					A.a[offset + 3] = -4.0 * robin1;
					A.a[offset + 4] = 1.0 + 25.0 * robin1 / 12.0;
					A.ja[offset] = BASE + IDX(i, NzInterior - 3);
					A.ja[offset + 1] = BASE + IDX(i, NzInterior - 2);
					A.ja[offset + 2] = BASE + IDX(i, NzInterior - 1);
					A.ja[offset + 3] = BASE + IDX(i, NzInterior);
					A.ja[offset + 4] = BASE + IDX(i, NzInterior + 1);
					break;
				}
				f[IDX(i, j)] = uInf;
			}
		}

		// At this point we have filled:
		offset = (2 + 2) + 2 * NzInterior
			+ 2 + 7 + 8 * (NzInterior - 2) + 9 + n_robin
			+ (NrInterior - 2) * (2 + 8 + 9 * (NzInterior - 2) + 10 + n_robin);

		// offset = (2 + 2) + 2 * NzInterior 
		// + 2 + 7 + 8 * (NzInterior - 2) + 9 + n_robin
		// + (NrInterior - 2) * (2 + 8 + 9 * (NzInterior - 2) + 10 + n_robin).

		// Interior bottom-right corner with equatorial symmetry.
		i = NrInterior;
		r = (double)i - 0.5;
		ir = 1.0 / r;

		j = 0;
		A.ia[IDX(NrInterior, 0)] = BASE + offset;
		A.a[offset] = 1.0;
		A.a[offset + 1] = -(double)z_sym;
		A.ja[offset] = BASE + IDX(NrInterior, 0);
		A.ja[offset + 1] = BASE + IDX(NrInterior, 1);
		f[IDX(NrInterior, 0)] = 0.0;
		offset += 2;

		// offset = (2 + 2) + 2 * NzInterior 
		// + 2 + 7 + 8 * (NzInterior - 2) + 9 + n_robin
		// + (NrInterior - 2) * (2 + 8 + 9 * (NzInterior - 2) + 10 + n_robin)
		// + 2.

		// Semi-onesided and centered difference.
		//
		//             o
		//             |
		//             o
		//             |
		// o--o--o--o--x--o
		//             |
		//             o
		//
		// Recall:
		// - Semi-onesided on r:
		// u''(x) = (10u(i+1)-15u(i)-4u(i-1)+14u(i-2)-6u(i-3)+u(i-4)/12h**2
		// u' (x) = (3u(i+1)+10u(i)-18u(i-1)+6u(i-2)-u(i-3))/12h
		// - Semi-onesided on z:
		// u''(x) = (-u(i+3)+4u(i+2)+6u(i+1)-20u(i)-11u(i-1))/12h**2.
		// 
		// Thus:
		// 
		// u(i-4, j  ) : 1/12
		// u(i-3, j  ) : -6/12 - 1h/12r
		// u(i-2, j  ) : 14/12 + 6h/12r
		// u(i-1, j  ) : -4/12 - 18h/12r
		// u(i  , j-1) : 16/12
		// u(i  , j  ) : s(i,j) * h**2 -15/12 - 30/12 +10h/12r
		// u(i  , j+1) : 16/12 + (z_sym)*(-1/12)
		// u(i  , j+2) : -1/12
		// u(i+1, j  ) : 10/12 + 3h/12r
		//
		j = 1;
		A.ia[IDX(i, j)] = BASE + offset;
		A.a[offset] = twelfth * zor;
		A.a[offset + 1] = (-0.5 - twelfth * ir) * zor;
		A.a[offset + 2] = (7.0 * sixth + 0.5 * ir) * zor;
		A.a[offset + 3] = (-third - 1.5 * ir) * zor;
		A.a[offset + 4] = 4.0 * third * roz;
		A.a[offset + 5] = dr * dz * s[IDX(i, j)] - 2.5 * roz - 1.25 * zor + 5.0 * sixth * ir * zor;
		A.a[offset + 6] = (4.0 * third + (double)z_sym * (-twelfth)) * roz;
		A.a[offset + 7] = -twelfth * roz;
		A.a[offset + 8] = (5.0 * sixth + 0.25 * ir) * zor;
		A.ja[offset] = BASE + IDX(i - 4, j);
		A.ja[offset + 1] = BASE + IDX(i - 3, j);
		A.ja[offset + 2] = BASE + IDX(i - 2, j);
		A.ja[offset + 3] = BASE + IDX(i - 1, j);
		A.ja[offset + 4] = BASE + IDX(i, j - 1);
		A.ja[offset + 5] = BASE + IDX(i, j);
		A.ja[offset + 6] = BASE + IDX(i, j + 1);
		A.ja[offset + 7] = BASE + IDX(i, j + 2);
		A.ja[offset + 8] = BASE + IDX(i + 1, j);
		f[IDX(i, j)] *= dr * dz;
		offset += 9;

		// offset = (2 + 2) + 2 * NzInterior 
		// + 2 + 7 + 8 * (NzInterior - 2) + 9 + n_robin
		// + (NrInterior - 2) * (2 + 8 + 9 * (NzInterior - 2) + 10 + n_robin).
		// + 2 + 9.

		// Set temporary offset.
		t_offset = offset;

		#pragma omp parallel shared(A, f) private(offset)
		{
			#pragma omp for schedule(guided)
			for (j = 2; j < NzInterior; j++)
			{
				// Eac iteration fills 10 elements.
				offset = t_offset + (j - 2) * 10;
				// Semi-onesided and centered difference.
				//
				//             o
				//             |
				//             o
				//             |
				// o--o--o--o--x--o
				//             |
				//             o
				//             |
				//             o
				//
				// Recall:
				// - Semi-onesided on r:
				// u''(x) = (10u(i+1)-15u(i)-4u(i-1)+14u(i-2)-6u(i-3)+u(i-4)/12h**2
				// u' (x) = (3u(i+1)+10u(i)-18u(i-1)+6u(i-2)-u(i-3))/12h
				// - Centered on z:
				// u''(x) = (-(u(i+2)+u(i-2))+16(u(i+1)+u(i-1))-30u(i))/12h**2
				// 
				// Thus:
				// 
				// u(i-4, j  ) : 1/12
				// u(i-3, j  ) : -6/12 - 1h/12r
				// u(i-2, j  ) : 14/12 + 6h/12r
				// u(i-1, j  ) : -4/12 - 18h/12r
				// u(i  , j-2) : -1/12
				// u(i  , j-1) : 16/12
				// u(i  , j  ) : s(i,j) * h**2 -15/12 + 10h/12r - 30/12
				// u(i  , j+1) : 16/12
				// u(i  , j+2) : -1/12
				// u(i+1, j  ) : 10/12 + 3h/12r
				//
				A.ia[IDX(i, j)] = BASE + offset;
				A.a[offset] = twelfth * zor;
				A.a[offset + 1] = (-0.5 - twelfth * ir) * zor;
				A.a[offset + 2] = (7.0 * sixth + 0.5 * ir) * zor;
				A.a[offset + 3] = (-third - 1.5 * ir) * zor;
				A.a[offset + 4] = -twelfth * roz;
				A.a[offset + 5] = 4.0 * third * roz;
				A.a[offset + 6] = dr * dz * s[IDX(i, j)] - 2.5 * roz - 1.25 * zor + 5.0 * sixth * ir * zor;
				A.a[offset + 7] = 4.0 * third * roz;
				A.a[offset + 8] = -twelfth * roz;
				A.a[offset + 9] = (5.0 * sixth + 0.25 * ir) * zor;
				A.ja[offset] = BASE + IDX(i - 4, j);
				A.ja[offset + 1] = BASE + IDX(i - 3, j);
				A.ja[offset + 2] = BASE + IDX(i - 2, j);
				A.ja[offset + 3] = BASE + IDX(i - 1, j);
				A.ja[offset + 4] = BASE + IDX(i, j - 2);
				A.ja[offset + 5] = BASE + IDX(i, j - 1);
				A.ja[offset + 6] = BASE + IDX(i, j);
				A.ja[offset + 7] = BASE + IDX(i, j + 1);
				A.ja[offset + 8] = BASE + IDX(i, j + 2);
				A.ja[offset + 9] = BASE + IDX(i + 1, j);
				f[IDX(i, j)] *= dr * dz;
			}
		}

		// At this point we have filled:
		offset = (2 + 2) + 2 * NzInterior
			+ 2 + 7 + 8 * (NzInterior - 2) + 9 + n_robin
			+ (NrInterior - 2) * (2 + 8 + 9 * (NzInterior - 2) + 10 + n_robin)
			+ 2 + 9 + 10 * (NzInterior - 2);

		// offset = (2 + 2) + 2 * NzInterior 
		// + 2 + 7 + 8 * (NzInterior - 2) + 9 + n_robin
		// + (NrInterior - 2) * (2 + 8 + 9 * (NzInterior - 2) + 10 + n_robin).
		// + 2 + 9 + 10 * (NzInterior - 2)

		// Fill interior top-right corner.
		//
		//             o
		//             |
		// o--o--o--o--x--o
		//             |
		//             o
		//             |
		//             o
		//             |
		//             o
		//             |
		//             o
		//
		// Recall:
		// - Semi-onesided on r:
		// u''(x) = (10u(i+1)-15u(i)-4u(i-1)+14u(i-2)-6u(i-3)+u(i-4)/12h**2
		// u' (x) = (3u(i+1)+10u(i)-18u(i-1)+6u(i-2)-u(i-3))/12h
		// - Semi-onesided on z:
		// u''(x) = (10u(i+1)-15u(i)-4u(i-1)+14u(i-2)-6u(i-3)+u(i-4)/12h**2
		//
		// Thus:
		//
		// u(i-4, j  ) : 1/12
		// u(i-3, j  ) : -6/12 - 1h/12r
		// u(i-2, j  ) : 14/12 + 6h/12r
		// u(i-1, j  ) : -4/12 - 18h/12r
		// u(i  , j-4) : 1/12
		// u(i  , j-3) : -6/12
		// u(i  , j-2) : 14/12
		// u(i  , j-1) : -4/12
		// u(i  , j  ) : s(i,j) * h**2 - 30/12 + 10h/12r
		// u(i  , j+1) : 10/12
		// u(i+1, j  ) : 10/12 + 3h/12r
		//
		j = NzInterior;
		A.ia[IDX(i, j)] = BASE + offset;
		A.a[offset] = twelfth * zor;
		A.a[offset + 1] = (-0.5 - twelfth * ir) * zor;
		A.a[offset + 2] = (7.0 * sixth + 0.5 * ir) * zor;
		A.a[offset + 3] = (-third - 1.5 * ir) * zor;
		A.a[offset + 4] = twelfth * roz;
		A.a[offset + 5] = -0.5 * roz;
		A.a[offset + 6] = 7.0 * sixth * roz;
		A.a[offset + 7] = -third * roz;
		A.a[offset + 8] = dr * dz * s[IDX(i, j)] - 1.25 * (roz + zor) + 5.0 * sixth * ir * zor;
		A.a[offset + 9] = 5.0 * sixth * roz;
		A.a[offset + 10] = (5.0 * sixth + 0.25 * ir) * zor;
		A.ja[offset] = BASE + IDX(i - 4, j);
		A.ja[offset + 1] = BASE + IDX(i - 3, j);
		A.ja[offset + 2] = BASE + IDX(i - 2, j);
		A.ja[offset + 3] = BASE + IDX(i - 1, j);
		A.ja[offset + 4] = BASE + IDX(i, j - 4);
		A.ja[offset + 5] = BASE + IDX(i, j - 3);
		A.ja[offset + 6] = BASE + IDX(i, j - 2);
		A.ja[offset + 7] = BASE + IDX(i, j - 1);
		A.ja[offset + 8] = BASE + IDX(i, j);
		A.ja[offset + 9] = BASE + IDX(i, j + 1);
		A.ja[offset + 10] = BASE + IDX(i + 1, j);
		f[IDX(i, j)] *= dr * dz;
		offset += 11;

		// offset = (2 + 2) + 2 * NzInterior 
		// + 2 + 7 + 8 * (NzInterior - 2) + 9 + n_robin
		// + (NrInterior - 2) * (2 + 8 + 9 * (NzInterior - 2) + 10 + n_robin).
		// + 2 + 9 + 10 * (NzInterior - 2) + 11

		j = NzInterior + 1;
		z = (double)j - 0.5;
		rr2 = r * r * roz * roz + z * z;
		A.ia[IDX(i, j)] = BASE + offset;
		switch (robin)
		{
		case 1:
			robin1 = rr2 / z;
			A.a[offset] = robin1 / 4.0;
			A.a[offset + 1] = -4.0 * robin1 / 3.0;
			A.a[offset + 2] = 3.0 * robin1;
			A.a[offset + 3] = -4.0 * robin1;
			A.a[offset + 4] = 1.0 + 25.0 * robin1 / 12.0;
			A.ja[offset] = BASE + IDX(i, NzInterior - 3);
			A.ja[offset + 1] = BASE + IDX(i, NzInterior - 2);
			A.ja[offset + 2] = BASE + IDX(i, NzInterior - 1);
			A.ja[offset + 3] = BASE + IDX(i, NzInterior);
			A.ja[offset + 4] = BASE + IDX(i, NzInterior + 1);
			break;
		}
		f[IDX(i, j)] = uInf;
		offset += n_robin;

		// offset = (2 + 2) + 2 * NzInterior 
		// + 2 + 7 + 8 * (NzInterior - 2) + 9 + n_robin
		// + (NrInterior - 2) * (2 + 8 + 9 * (NzInterior - 2) + 10 + n_robin).
		// + 2 + 9 + 10 * (NzInterior - 2) + 11 + n_robin

		// Do bottom-right corner.
		i = NrInterior + 1;
		r = (double)i - 0.5;

		j = 0;
		A.ia[IDX(NrInterior + 1, 0)] = BASE + offset;
		A.a[offset] = 1.0;
		A.a[offset + 1] = -(double)z_sym;
		A.ja[offset] = BASE + IDX(NrInterior + 1, 0);
		A.ja[offset + 1] = BASE + IDX(NrInterior + 1, 1);
		f[IDX(NrInterior + 1, 0)] = 0.0;
		offset += 2;

		// offset = (2 + 2 + 2) + 2 * NzInterior 
		// + 2 + 7 + 8 * (NzInterior - 2) + 9 + n_robin
		// + (NrInterior - 2) * (2 + 8 + 9 * (NzInterior - 2) + 10 + n_robin).
		// + 2 + 9 + 10 * (NzInterior - 2) + 11 + n_robin

		// Set temporary offset.
		t_offset = offset;

		#pragma omp parallel shared(A, f) private(offset, z, rr2,\
		robin1, robin2, robin3)
		{
			#pragma omp for schedule(guided)
			for (j = 1; j < NzInterior + 1; j++)
			{
				// Each iteration fills n_robin elements.
				offset = t_offset + n_robin * (j - 1);

				// Z coordinate.
				z = (double)j - 0.5;
				rr2 = r * r + z * z * zor * zor;
				A.ia[IDX(i, j)] = BASE + offset;
				switch (robin)
				{
				case 1:
					robin1 = rr2 / r;
					A.a[offset] = robin1 / 4.0;
					A.a[offset + 1] = -4.0 * robin1 / 3.0;
					A.a[offset + 2] = 3.0 * robin1;
					A.a[offset + 3] = -4.0 * robin1;
					A.a[offset + 4] = 1.0 + 25.0 * robin1 / 12.0;
					A.ja[offset] = BASE + IDX(NrInterior - 3, j);
					A.ja[offset + 1] = BASE + IDX(NrInterior - 2, j);
					A.ja[offset + 2] = BASE + IDX(NrInterior - 1, j);
					A.ja[offset + 3] = BASE + IDX(NrInterior, j);
					A.ja[offset + 4] = BASE + IDX(NrInterior + 1, j);
					break;
				}
				f[IDX(i, j)] = uInf;
			}
		}

		offset = (2 + 2 + 2) + 2 * NzInterior
			+ 2 + 7 + 8 * (NzInterior - 2) + 9 + n_robin
			+ (NrInterior - 2) * (2 + 8 + 9 * (NzInterior - 2) + 10 + n_robin)
			+ 2 + 9 + 10 * (NzInterior - 2) + 11 + n_robin
			+ n_robin * NzInterior;

		// offset = (2 + 2 + 2) + 2 * NzInterior 
		// + 2 + 7 + 8 * (NzInterior - 2) + 9 + n_robin
		// + (NrInterior - 2) * (2 + 8 + 9 * (NzInterior - 2) + 10 + n_robin).
		// + 2 + 9 + 10 * (NzInterior - 2) + 11 + n_robin
		// + n_robin * NzInterior

		// Finally, fill top corner with Robin
		j = NzInterior + 1;
		z = (double)j - 0.5;
		rrodrr = sqrt((r * r * dr * dr + z * z * dz * dz) / (dr * dr + dz * dz));
		A.ia[IDX(i, j)] = BASE + offset;
		switch (robin)
		{
		case 1:
			A.a[offset] = rrodrr / 4.0;
			A.a[offset + 1] = -4.0 * rrodrr / 3.0;
			A.a[offset + 2] = 3.0 * rrodrr;
			A.a[offset + 3] = -4.0 * rrodrr;
			A.a[offset + 4] = 1.0 + 25.0 * rrodrr / 12.0;
			A.ja[offset] = BASE + IDX(NrInterior - 3, NzInterior - 3);
			A.ja[offset + 1] = BASE + IDX(NrInterior - 2, NzInterior - 2);
			A.ja[offset + 2] = BASE + IDX(NrInterior - 1, NzInterior - 1);
			A.ja[offset + 3] = BASE + IDX(NrInterior, NzInterior);
			A.ja[offset + 4] = BASE + IDX(NrInterior + 1, NzInterior + 1);
			break;
		}
		f[IDX(i, j)] = uInf;
		offset += n_robin;

		// offset = (2 + 2 + 2 + n_robin) + 2 * NzInterior 
		// + 2 + 7 + 8 * (NzInterior - 2) + 9 + n_robin
		// + (NrInterior - 2) * (2 + 8 + 9 * (NzInterior - 2) + 10 + n_robin).
		// + 2 + 9 + 10 * (NzInterior - 2) + 11 + n_robin
		// + n_robin * NzInterior

		// Fill last element.
		A.ia[IDX(NrInterior + 1, NzInterior + 1) + 1] = BASE + nnz;
	}

#ifdef DEBUG
	csr_print(A, "A_a.asc", "A_ia.asc", "A_ja.asc");
#endif

	// All done.
	return;
}
