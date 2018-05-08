// Global header.
// One-based indexing BASE is defined in this header.
#include "tools.h"

// Print CSR matrix for debug.
#undef DEBUG

// Nonzero calculator.
int nnz_general_elliptic(const int NrInterior, const int NzInterior, const int order, const int robin)
{
	int nnz;

	// Number of elements per Robin discretization.
	int n_robin = robin + order;

	// Second order general elliptic equation.
	if (order == 2)
	{
		nnz = 9 * NrInterior * NzInterior
			+ (2 + n_robin) * (NrInterior + NzInterior)
			+ (2 + 2 + 2 + n_robin);
	}
	// Fourth order general elliptic equation.
	else if (order == 4)
	{
		nnz = 17 * (NrInterior - 2) * (NzInterior - 2)
			+ (16 + 26) * (NrInterior + NzInterior - 4)
			+ (14 + 21 + 21 + 27)
			+ (2 + n_robin) * (NrInterior + NzInterior)
			+ (2 + 2 + 2 + n_robin);
	}
	return nnz;
}

// Write CSR matrix for the general elliptic equation.
void csr_gen_general_elliptic(csr_matrix A,	// CSR matrix structure.
	const int NrInterior,			// Number of r interior points.	
	const int NzInterior,			// Number of z interior points.
	const int order,			// Finite difference order: 2 or 4.
	const double dr,			// Spatial step in r.
	const double dz,			// Spatial step in z.
	const double *ell_a,			// Coefficient of (d^2/dr^2)
	const double *ell_b,			// Coefficient of (d^2/drdz)
	const double *ell_c,			// Coefficient of (d^2/dz^)
	const double *ell_d,			// Coefficient of (d/dr)
	const double *ell_e,			// Coefficient of (d/dz)
	const double *ell_s,			// Linear source.
	double *ell_f,				// RHS.
	const double uInf,			// Value at infinity.
	const int robin,			// Robin BC type: 1, 2, 3.
	const int r_sym,			// R symmetry: 1(even), -1(odd).
	const int z_sym)			// Z symmetry: 1(even), -1(odd).
{
	// Constant numbers.
	const double third = 1.0 / 3.0;
	const double sixth = 1.0 / 6.0;
	const double twelfth = 1.0 / 12.0;

	// Grid extensions.
	int NrTotal = NrInterior + 2;
	int NzTotal = NzInterior + 2;

	// Number of nonzero elements.
	int nnz = nnz_general_elliptic(NrInterior, NzInterior, order, robin);

	// Number of Robin elements.
	int n_robin = robin + order;

	// Number of elements we have filled in.
	int offset = 0;

	// Auxiliary variables.
	double r, z, rr2, rrodrr;
	int i, j, t_offset;
	double aux_a, aux_b, aux_c, aux_d, aux_e, aux_s;
	double robin1, robin2, robin3;

	// Ratios for variable step size.
	double roz = dr / dz;
	double zor = dz / dr;

	// SECOND-ORDER GENERAL ELLIPTIC EQUATION.
	if (order == 2)
	{
        // Start offset at zero.
		offset = 0;

		// Lower-left corner: diagonal symmetry.
		A.ia[IDX(0, 0)] = BASE + offset;
		A.a[offset] = 1.0;
		A.a[offset + 1] = -(double)(r_sym * z_sym);
		A.ja[offset] = BASE + IDX(0, 0);
		A.ja[offset + 1] = BASE + IDX(1, 1);
		ell_f[IDX(0, 0)] = 0.0;
		offset += 2;

		// Set temporary offset.
		t_offset = offset;

		// Fill left-boundary using axis symmetry.
		#pragma omp parallel shared(A, ell_f) private(offset)
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
				ell_f[IDX(0, j)] = 0.0;
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
		ell_f[IDX(0, NzInterior + 1)] = 0.0;
		offset += 2;

		// Set temporary offset.
		t_offset = offset;

		// Now come the interior points plus the top and bottom boundaries with
		// Robin and equatorial symmetry respectively.
		#pragma omp parallel shared(A, ell_f) private(offset, j, r, z, rr2,\
		aux_a, aux_b, aux_c, aux_d, aux_e, aux_s,\
		robin1, robin2, robin3)
		{
			#pragma omp for schedule(guided)
			for (i = 1; i < NrInterior + 1; i++)
			{
				// Each iteration of i loop will fill 9 * NzInterior + (2 + n_robin) values.
				offset = t_offset + (i - 1) * (9 * NzInterior + 2 + n_robin);

				// R coordinate.
				r = (double)i - 0.5;

				// Do bottom boundary first with equatorial symmetry.
				A.ia[IDX(i, 0)] = BASE + offset;
				A.a[offset] = 1.0;
				A.a[offset + 1] = -(double)z_sym;
				A.ja[offset] = BASE + IDX(i, 0);
				A.ja[offset + 1] = BASE + IDX(i, 1);
				ell_f[IDX(i, 0)] = 0.0;
				offset += 2;

				// Now loop over interior points.
				for (j = 1; j < NzInterior + 1; j++)
				{
					// Fetch values, notice the dr * dz rescaling.
					aux_a = ell_a[IDX(i, j)] * zor;
					aux_b = 0.25 * ell_b[IDX(i, j)];
					aux_c = ell_c[IDX(i, j)] * roz;
					aux_d = 0.5 * dz * ell_d[IDX(i, j)];
					aux_e = 0.5 * dr * ell_e[IDX(i, j)];
					aux_s = dr * dz * ell_s[IDX(i, j)];
					// Row begins at offset.
					A.ia[IDX(i, j)] = BASE + offset;
					// Values.
					A.a[offset] = aux_b;
					A.a[offset + 1] = aux_a - aux_d;
					A.a[offset + 2] = -aux_b;
					A.a[offset + 3] = aux_c - aux_e;
					A.a[offset + 4] = aux_s - 2.0 * (aux_a + aux_c);
					A.a[offset + 5] = aux_c + aux_e;
					A.a[offset + 6] = -aux_b;
					A.a[offset + 7] = aux_a + aux_d;
					A.a[offset + 8] = aux_b;
					// Columns.
					A.ja[offset] = BASE + IDX(i - 1, j - 1);
					A.ja[offset + 1] = BASE + IDX(i - 1, j);
					A.ja[offset + 2] = BASE + IDX(i - 1, j + 1);
					A.ja[offset + 3] = BASE + IDX(i, j - 1);
					A.ja[offset + 4] = BASE + IDX(i, j);
					A.ja[offset + 5] = BASE + IDX(i, j + 1);
					A.ja[offset + 6] = BASE + IDX(i + 1, j - 1);
					A.ja[offset + 7] = BASE + IDX(i + 1, j);
					A.ja[offset + 8] = BASE + IDX(i + 1, j + 1);
					// Multiply RHS by dr * dz.
					ell_f[IDX(i, j)] *= dr * dz;
					// Increase offset.
					offset += 9;
				}

				// Now fill the top boundary with Robin.
				j = NzInterior + 1;
				// Z coordinate.
				z = (double)j - 0.5;
				// Radial coordinate: overal division by dz**2.
				rr2 = r * r * roz * roz + z * z;
				A.ia[IDX(i, NzInterior + 1)] = BASE + offset;
				switch (robin)
				{
					case 1:
						robin1 = rr2 / z;
						A.a[offset] = 0.5 * robin1;
						A.a[offset + 1] = -2.0 * robin1;
						A.a[offset + 2] = 1.0 + 1.5 * robin1;
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
				ell_f[IDX(i, NzInterior + 1)] = uInf;
			}
		}

		// At this point we have now filled:
		offset = 4 + 2 * NzInterior + 9 * NrInterior * NzInterior + (2 + n_robin) * NrInterior;

		// Lower-right corner: equatorial symmetry.
		A.ia[IDX(NrInterior + 1, 0)] = BASE + offset;
		A.a[offset] = 1.0;
		A.a[offset + 1] = -(double)z_sym;
		A.ja[offset] = BASE + IDX(NrInterior + 1, 0);
		A.ja[offset + 1] = BASE + IDX(NrInterior + 1, 1);
		ell_f[IDX(NrInterior + 1, 0)] = 0.0;
		offset += 2;

		// Set temporary offset.
		t_offset = offset;

		// Robin boundary.
		r = (double)NrInterior + 0.5;
		#pragma omp parallel shared(A, ell_f) private(offset, z, rr2,\
		robin1, robin2, robin3)
		{
			#pragma omp for schedule(guided)
			for (j = 1; j < NzInterior + 1; j++)
			{
                		// Overall division by dr**2.
				z = (double)j - 0.5;
				rr2 = r * r + z * z * zor * zor;

				// Each iteration of the loop fills n_robin elements.
				offset = t_offset + n_robin * (j - 1);

				A.ia[IDX(NrInterior + 1, j)] = BASE + offset;
				switch (robin)
				{
					case 1:
						robin1 = (rr2 / r);
						A.a[offset] = 0.5 * robin1;
						A.a[offset + 1] = -2.0 * robin1;
						A.a[offset + 2] = 1.0 + 1.5 * robin1;
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
				ell_f[IDX(NrInterior + 1, j)] = uInf;
			}
		}

		// At this point, we have now filled:
		offset = 6 + (2 + n_robin) * (NzInterior + NrInterior) + 9 * NrInterior * NzInterior;

		// Upper-right corner: fill with Robin.
		r = (double)NrInterior + 0.5;
		z = (double)NzInterior + 0.5;
        	// Radial coordinate over radial step.
		rrodrr = sqrt((r * r * dr * dr + z * z * dz * dz) / (dr * dr + dz * dz));
		A.ia[IDX(NrInterior + 1, NzInterior + 1)] = BASE + offset;
		switch (robin)
		{
			case 1:
				robin1 = rrodrr;
				A.a[offset] = 0.5 * robin1;
				A.a[offset + 1] = -2.0 * robin1;
				A.a[offset + 2] = 1.0 + 1.5 * robin1;
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
		ell_f[IDX(NrInterior + 1, NzInterior + 1)] = uInf;

		// Assert fill-in.
		assert(nnz == offset);

		// Fill last element of row offsets.
		A.ia[IDX(NrInterior + 1, NzInterior + 1) + 1] = BASE + nnz;
	}
	// FOURTH-ORDER GENERAL ELLIPTIC EQUATION
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
		ell_f[IDX(0, 0)] = 0.0;
		offset += 2;

		// offset = 2.

		// Set temporary offset.
		t_offset = offset;

		// Fill left-boundary using axis symmetry.
		#pragma omp parallel shared(A, ell_f) private(offset)
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
				ell_f[IDX(0, j)] = 0.0;
			}
		}

#ifdef DEBUG
		printf("CSR FILL-IN: Stage 1.\n");
#endif

		// We have now filled:
		offset = 2 + 2 * NzInterior;

		// offset = 2 + NzInterior * 2.

		// Upper-left corner: also axis symmetry.
		A.ia[IDX(0, NzInterior + 1)] = BASE + offset;
		A.a[offset] = 1.0;
		A.a[offset + 1] = -(double)r_sym;
		A.ja[offset] = BASE + IDX(0, NzInterior + 1);
		A.ja[offset + 1] = BASE + IDX(1, NzInterior + 1);
		ell_f[IDX(0, NzInterior + 1)] = 0.0;
		offset += 2;

		// offset = (2 + 2) + 2 * NzInterior.

		// Left-interior points with semi-one sided finite difference.
		i = 1;
		r = 0.5;

		// Bottom boundary with equatorial symmetry.
		j = 0;
		A.ia[IDX(1, 0)] = BASE + offset;
		A.a[offset] = 1.0;
		A.a[offset + 1] = -(double)z_sym;
		A.ja[offset] = BASE + IDX(1, 0);
		A.ja[offset + 1] = BASE + IDX(1, 1);
		ell_f[IDX(1, 0)] = 0.0;
		offset += 2;

		// offset = (2 + 2) + 2 * NzInterior + 2.

		// Stencil for lower-left corner.
		//
		//     o   o   o
		//     |   | / 
		// o   o   o---o
		//   \ | / 
		// o---x---o---o
		//   / | \
		// o   o   o
		//
		// Thus:
		// 
		// u(i-1, j-1) : b/3
		// u(i-1, j  ) : 2(2a - d)/3
		// u(i-1, j+1) : -b/3
		// u(i  , j-1) : 2(2c - e)/3
		// u(i  , j  ) : -5(a + c)/2 + s
		// u(i  , j+1) : 2(2c + e)/3 + z_sym*((-c + e)/12)
		// u(i  , j+2) : -(c + e)/12
		// u(i+1, j-1) : -b/3
		// u(i+1, j  ) : 2(2a + d)/3 + r_sym*((-a + d)/12)
		// u(i+1, j+1) : b/3 + r_sym*z_sym*(-b/48)
		// u(i+1, j+2) : r_sym*(b/48)
		// u(i+2, j  ) : -(a + d)/12
		// u(i+2, j+1) : z_sym*(b/48)
		// u(i+2, j+2) : -b/48
		//
		j = 1;
		aux_a = ell_a[IDX(i, j)] * zor;
		aux_b = ell_b[IDX(i, j)];
		aux_c = ell_c[IDX(i, j)] * roz;
		aux_d = dz * ell_d[IDX(i, j)];
		aux_e = dr * ell_e[IDX(i, j)];
		aux_s = dr * dz * ell_s[IDX(i, j)];

		A.ia[IDX(1, 1)] = BASE + offset;
		A.a[offset] = aux_b * third;
		A.a[offset + 1] = 2.0 * (2.0 * aux_a - aux_d) * third;
		A.a[offset + 2] = -aux_b * third;
		A.a[offset + 3] = 2.0 * (2.0 * aux_c - aux_e) * third;
		A.a[offset + 4] = -2.5 * (aux_a + aux_c) + aux_s;
		A.a[offset + 5] = 2.0 * (2.0 * aux_c + aux_e) * third + (double)z_sym * ((-aux_c + aux_e) * twelfth);
		A.a[offset + 6] = -(aux_c + aux_e) * twelfth;
		A.a[offset + 7] = -aux_b * third;
		A.a[offset + 8] = 2.0 * (2.0 * aux_a + aux_d) * third + (double)r_sym * ((-aux_a + aux_d) * twelfth);
		A.a[offset + 9] = aux_b * third + (double)(r_sym * z_sym) * (-aux_b / 48.0);
		A.a[offset + 10] = (double)r_sym * (aux_b / 48.0);
		A.a[offset + 11] = -(aux_a + aux_d) * twelfth;
		A.a[offset + 12] = (double)z_sym * (aux_b / 48.0);
		A.a[offset + 13] = -aux_b / 48.0;
		A.ja[offset] = BASE + IDX(0, 0);
		A.ja[offset + 1] = BASE + IDX(0, 1);
		A.ja[offset + 2] = BASE + IDX(0, 2);
		A.ja[offset + 3] = BASE + IDX(1, 0);
		A.ja[offset + 4] = BASE + IDX(1, 1);
		A.ja[offset + 5] = BASE + IDX(1, 2);
		A.ja[offset + 6] = BASE + IDX(1, 3);
		A.ja[offset + 7] = BASE + IDX(2, 0);
		A.ja[offset + 8] = BASE + IDX(2, 1);
		A.ja[offset + 9] = BASE + IDX(2, 2);
		A.ja[offset + 10] = BASE + IDX(2, 3);
		A.ja[offset + 11] = BASE + IDX(3, 1);
		A.ja[offset + 12] = BASE + IDX(3, 2);
		A.ja[offset + 13] = BASE + IDX(3, 3);
		ell_f[IDX(1, 1)] *= dr * dz;
		offset += 14;

		// offset = (2 + 2) +  2 * NzInterior + (2 + 14).

		// Set temporary offset.
		t_offset = offset;

		// Stencil for left strip.
		//
		//     o   o   o
		//     | /   / 
		// o   o   o    
		//   \ | / 
		// o---x---o---o
		//   / | \
		// o   o   o
		//     | \   \
		//     o   o   o
		//
		// Thus:
		// 
		// u(i-1, j-1) : b/3
		// u(i-1, j  ) : 2(2a - d)/3
		// u(i-1, j+1) : -b/3
		// u(i  , j-2) : (-c + e)/12
		// u(i  , j-1) : 2(2c - e)/3
		// u(i  , j  ) : -5(a + c)/2 + s
		// u(i  , j+1) : 2(2c + e)/3
		// u(i  , j+2) : -(c + e)/12
		// u(i+1, j-2) : r_sym*(-b/48)
		// u(i+1, j-1) : -b/3
		// u(i+1, j  ) : 2(2a + d)/3 + r_sym*((-a + d)/12)
		// u(i+1, j+1) : b/3
		// u(i+1, j+2) : r_sym*(b/48)
		// u(i+2, j-2) : b/48
		// u(i+2, j  ) : -(a + d)/12
		// u(i+2, j+2) : -b/48
		// 
		#pragma omp parallel shared(A, ell_f) private(offset,\
		aux_a, aux_b, aux_c, aux_d, aux_e, aux_s)
		{
			#pragma omp for schedule(guided)
			for (j = 2; j < NzInterior; j++)
			{
				// Each iterations fills 16 elements.
				offset = t_offset + 16 * (j - 2);
				aux_a = ell_a[IDX(i, j)] * zor;
				aux_b = ell_b[IDX(i, j)];
				aux_c = ell_c[IDX(i, j)] * roz;
				aux_d = dz * ell_d[IDX(i, j)];
				aux_e = dr * ell_e[IDX(i, j)];
				aux_s = dr * dz * ell_s[IDX(i, j)];

				A.ia[IDX(1, j)] = BASE + offset;
				A.a[offset] = aux_b * third;
				A.a[offset + 1] = 2.0 * (2.0 * aux_a - aux_d) * third;
				A.a[offset + 2] = -aux_b * third;
				A.a[offset + 3] = (-aux_c + aux_e) * twelfth;
				A.a[offset + 4] = 2.0 * (2.0 * aux_c - aux_e) * third;
				A.a[offset + 5] = -2.5 * (aux_a + aux_c) + aux_s;
				A.a[offset + 6] = 2.0 * (2.0 * aux_c + aux_e) * third;
				A.a[offset + 7] = -(aux_c + aux_e) * twelfth;
				A.a[offset + 8] = (double)r_sym * (-aux_b / 48.0);
				A.a[offset + 9] = -aux_b * third;
				A.a[offset + 10] = 2.0 * (2.0 * aux_a + aux_d) * third + (double)r_sym * ((-aux_a + aux_d) * twelfth);
				A.a[offset + 11] = aux_b * third;
				A.a[offset + 12] = (double)r_sym * (aux_b / 48.0);
				A.a[offset + 13] = aux_b / 48.0;
				A.a[offset + 14] = -(aux_a + aux_d) * twelfth;
				A.a[offset + 15] = -aux_b / 48.0;
				A.ja[offset] = BASE + IDX(i - 1, j - 1);
				A.ja[offset + 1] = BASE + IDX(i - 1, j);
				A.ja[offset + 2] = BASE + IDX(i - 1, j + 1);
				A.ja[offset + 3] = BASE + IDX(i, j - 2);
				A.ja[offset + 4] = BASE + IDX(i, j - 1);
				A.ja[offset + 5] = BASE + IDX(i, j);
				A.ja[offset + 6] = BASE + IDX(i, j + 1);
				A.ja[offset + 7] = BASE + IDX(i, j + 2);
				A.ja[offset + 8] = BASE + IDX(i + 1, j - 2);
				A.ja[offset + 9] = BASE + IDX(i + 1, j - 1);
				A.ja[offset + 10] = BASE + IDX(i + 1, j);
				A.ja[offset + 11] = BASE + IDX(i + 1, j + 1);
				A.ja[offset + 12] = BASE + IDX(i + 1, j + 2);
				A.ja[offset + 13] = BASE + IDX(i + 2, j - 2);
				A.ja[offset + 14] = BASE + IDX(i + 2, j);
				A.ja[offset + 15] = BASE + IDX(i + 2, j + 2);
				ell_f[IDX(1, j)] *= dr * dz;
			}
		}

#ifdef DEBUG
		printf("CSR FILL-IN: Stage 2.\n");
#endif

		// At this point we have filled:
		offset = (2 + 2) + 2 * NzInterior + 2 + 14 + 16 * (NzInterior - 2);

		// offset = (2 + 2) + 2 * NzInterior + 2 + 14 + 16 * (NzInterior - 2).


		// Stencil for top left corner.
		//
		// o   o   o   o
		//   \ | /   /
		// o---x---o---o
		//   / | \   \
		// o   o   o   o
		//   / | \   \
		// o   o   o   o
		//   / | \   \
		// o   o   o   o
		//     |      
		//     o        
		//
		// Thus:
		// 
		// u(i-1, j-3) : b/18
		// u(i-1, j-2) : -b/3
		// u(i-1, j-1) : b
		// u(i-1, j  ) : (12a - 5b - 6d)/9
		// u(i-1, j+1) : -b/6
		// u(i  , j-4) : c/12
		// u(i  , j-3) : -(6c + e)/12
		// u(i  , j-2) : (7c + 3e)/6
		// u(i  , j-1) : -(2c + 9e)/6
		// u(i  , j  ) : -5a/2 - 5c/4 + 5e/6 + s
		// u(i  , j+1) : 5c/6 + e/4
		// u(i+1, j-3) : -b/18 + r_sym*(-b/144)
		// u(i+1, j-2) : b/3 + r_sym*(b/24)
		// u(i+1, j-1) : -b + r_sym*(-b/8)
		// u(i+1, j  ) : (12a + 5b + 6d)/9 + r_sym*(-6a + 5b + 6d)/72
		// u(i+1, j+1) : b/6 + r_sym*(b/48)
		// u(i+2, j-3) : b/144
		// u(i+2, j-2) : -b/24
		// u(i+2, j-1) : b/8
		// u(i+2, j  ) : -(6a + 5b + 6d)/72
		// u(i+2, j+1) : -b/48
		//
		j = NzInterior;
		aux_a = ell_a[IDX(i, j)] * zor;
		aux_b = ell_b[IDX(i, j)];
		aux_c = ell_c[IDX(i, j)] * roz;
		aux_d = dz * ell_d[IDX(i, j)];
		aux_e = dr * ell_e[IDX(i, j)];
		aux_s = dr * dz * ell_s[IDX(i, j)];

		A.ia[IDX(1, j)] = BASE + offset;
		A.a[offset] = aux_b / 18.0;
		A.a[offset + 1] = -aux_b * third;
		A.a[offset + 2] = aux_b;
		A.a[offset + 3] = (12.0 * aux_a - 5.0 * aux_b - 6.0 * aux_d) / 9.0;
		A.a[offset + 4] = -aux_b * sixth;
		A.a[offset + 5] = aux_c * twelfth;
		A.a[offset + 6] = -(6.0 * aux_c + aux_e) * twelfth;
		A.a[offset + 7] = (7.0 * aux_c + 3.0 * aux_e) * sixth;
		A.a[offset + 8] = -(2.0 * aux_c + 9.0 * aux_e) * sixth;
		A.a[offset + 9] = -2.5 * aux_a - 1.25 * aux_c + 5.0 * sixth * aux_e + aux_s;
		A.a[offset + 10] = 5.0 * aux_c * sixth + aux_e * 0.25;
		A.a[offset + 11] = -aux_b / 18.0 + (double)r_sym * (-aux_b / 144.0);
		A.a[offset + 12] = aux_b * third + (double)r_sym * (aux_b / 24.0);
		A.a[offset + 13] = -aux_b + (double)r_sym * (-aux_b * 0.125);
		A.a[offset + 14] = (12.0 * aux_a + 5.0 * aux_b + 6.0 * aux_d) / 9.0 + (double)r_sym * ((-6.0 * aux_a + 5.0 * aux_b + 6.0 * aux_d) / 72.0);
		A.a[offset + 15] = aux_b * sixth + (double)r_sym * (aux_b / 48.0);
		A.a[offset + 16] = aux_b / 144.0;
		A.a[offset + 17] = -aux_b / 24.0;
		A.a[offset + 18] = aux_b * 0.125;
		A.a[offset + 19] = -(6.0 * aux_a + 5.0 * aux_b + 6.0 * aux_d) / 72.0;
		A.a[offset + 20] = -aux_b / 48.0;
		A.ja[offset] = BASE + IDX(i - 1, j - 3);
		A.ja[offset + 1] = BASE + IDX(i - 1, j - 2);
		A.ja[offset + 2] = BASE + IDX(i - 1, j - 1);
		A.ja[offset + 3] = BASE + IDX(i - 1, j);
		A.ja[offset + 4] = BASE + IDX(i - 1, j + 1);
		A.ja[offset + 5] = BASE + IDX(i, j - 4);
		A.ja[offset + 6] = BASE + IDX(i, j - 3);
		A.ja[offset + 7] = BASE + IDX(i, j - 2);
		A.ja[offset + 8] = BASE + IDX(i, j - 1);
		A.ja[offset + 9] = BASE + IDX(i, j);
		A.ja[offset + 10] = BASE + IDX(i, j + 1);
		A.ja[offset + 11] = BASE + IDX(i + 1, j - 3);
		A.ja[offset + 12] = BASE + IDX(i + 1, j - 2);
		A.ja[offset + 13] = BASE + IDX(i + 1, j - 1);
		A.ja[offset + 14] = BASE + IDX(i + 1, j);
		A.ja[offset + 15] = BASE + IDX(i + 1, j + 1);
		A.ja[offset + 16] = BASE + IDX(i + 2, j - 3);
		A.ja[offset + 17] = BASE + IDX(i + 2, j - 2);
		A.ja[offset + 18] = BASE + IDX(i + 2, j - 1);
		A.ja[offset + 19] = BASE + IDX(i + 2, j);
		A.ja[offset + 20] = BASE + IDX(i + 2, j + 1);
		ell_f[IDX(1, j)] *= dr * dz;
		offset += 21;

		// offset = (2 + 2) + 2 * NzInterior + 2 + 14 + 16 * (NzInterior - 2) + 21.
		j = NzInterior + 1;
		z = (double)j - 0.5;
        	// Overall division by dz**2.
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
			case 2:
				robin2 = (rr2 / z) * (rr2 / z);
				robin1 = (rr2 / z) * (4.0 - (r * roz / z) * (r * roz / z));
				A.a[offset] = -5.0 * robin2 / 12.0;
				A.a[offset + 1] = 61.0 * robin2 / 24.0 + 0.125 * robin1;
				A.a[offset + 2] = -(6.5 * robin2 + 2.0 * robin1 / 3.0);
				A.a[offset + 3] = 107.0 * robin2 / 12.0 + 1.5 * robin1;
				A.a[offset + 4] = -(77.0 * robin2 / 12.0 + 2.0 * robin1);
				A.a[offset + 5] = 1.0 + 1.875 * robin2 + 25.0 * robin1 / 24.0;
				A.ja[offset] = BASE + IDX(i, NzInterior - 4);
				A.ja[offset + 1] = BASE + IDX(i, NzInterior - 3);
				A.ja[offset + 2] = BASE + IDX(i, NzInterior - 2);
				A.ja[offset + 3] = BASE + IDX(i, NzInterior - 1);
				A.ja[offset + 4] = BASE + IDX(i, NzInterior);
				A.ja[offset + 5] = BASE + IDX(i, NzInterior + 1);
				break;
			case 3:
				robin3 = (rr2 / z) * (rr2 / z) * (rr2 / z);
				robin2 = (rr2 / z) * (rr2 / z) * (9.0 - 3.0 * (r * roz / z) * (r * roz / z));
				robin1 = (rr2 / z) * (18.0 + (r * roz / z) * (r * roz / z) * (-9.0 + 3.0 * (rr2 / (z * z))));
				A.a[offset] = 0.3125 * robin3;
				A.a[offset + 1] = -(13.0 * robin3 / 6.0 + 5.0 * robin2 / 36.0);
				A.a[offset + 2] = 307.0 * robin3 / 48.0 + 61.0 * robin2 / 72.0 + robin1 / 24.0;
				A.a[offset + 3] = -(31.0 * robin3 / 3.0 + 13.0 * robin2 / 6.0 + 2.0 * robin1 / 9.0);
				A.a[offset + 4] = 461.0 * robin3 / 48.0 + 107.0 * robin2 / 36.0 + 0.5 * robin1;
				A.a[offset + 5] = -(29.0 * robin3 / 6.0 + 77.0 * robin2 / 36.0 + 2.0 * robin1 / 3.0);
				A.a[offset + 6] = 1.0 + 49.0 * robin3 / 48.0 + 0.625 * robin2 + 25.0 * robin1 / 72.0;
				A.ja[offset] = BASE + IDX(i, NzInterior - 5);
				A.ja[offset + 1] = BASE + IDX(i, NzInterior - 4);
				A.ja[offset + 2] = BASE + IDX(i, NzInterior - 3);
				A.ja[offset + 3] = BASE + IDX(i, NzInterior - 2);
				A.ja[offset + 4] = BASE + IDX(i, NzInterior - 1);
				A.ja[offset + 5] = BASE + IDX(i, NzInterior);
				A.ja[offset + 6] = BASE + IDX(i, NzInterior + 1);
				break;
		}
		ell_f[IDX(1, NzInterior + 1)] = uInf;
		offset += n_robin;

		// offset = (2 + 2) + 2 * NzInterior 
		// + 2 + 14 + 16 * (NzInterior - 2) + 21 + n_robin.

		// Set temporary offset.
		t_offset = offset;

		#pragma omp parallel shared(A, ell_f) private(j, offset, r, z, rr2,\
		aux_a, aux_b, aux_c, aux_d, aux_e, aux_s,\
		robin1, robin2, robin3)
		{
			#pragma omp for schedule(guided)
			for (i = 2; i < NrInterior; i++)
			{
				// Each iteration fills 2 + 16 + 17 * (NzInterior - 2) + 26 + n_robin points.
				offset = t_offset + (i - 2) * (17 * NzInterior + 10 + n_robin);

				// R coordinate.
				r = (double)i - 0.5;

				// Equatorial symmetry.
				j = 0;
				A.ia[IDX(i, 0)] = BASE + offset;
				A.a[offset] = 1.0;
				A.a[offset + 1] = -(double)z_sym;
				A.ja[offset] = BASE + IDX(i, 0);
				A.ja[offset + 1] = BASE + IDX(i, 1);
				ell_f[IDX(i, 0)] = 0.0;
				offset += 2;

				// Stencil for lower strip.
				//
				// o       o       o
				//   \     |     /
				// o   o   o   o   o
				//   \   \ | /   /
				// o---o---x---o---o
				//       / | \
				//     o   o   o
				//
				// Thus:
				// 
				// u(i-2, j  ) : (-a + d)/12
				// u(i-2, j+1) : z_sym*(-b/48)
				// u(i-2, j+2) : b/48
				// u(i-1, j-1) : b/3
				// u(i-1, j  ) : 2(2a - d)/3
				// u(i-1, j+1) : -b/3
				// u(i  , j-1) : 2(2c - e)/3
				// u(i  , j  ) : -5(a + c)/2 + s
				// u(i  , j+1) : 2(2c + e)/3 + z_sym*((-c + e)/12)
				// u(i  , j+2) : -(c + e)/12
				// u(i+1, j-1) : -b/3
				// u(i+1, j  ) : 2(2a + d)/3
				// u(i+1, j+1) : b/3
				// u(i+2, j  ) : -(a + d)/12
				// u(i+2, j+1) : z_sym*(b/48)
				// u(i+2, j+2) : -b/48
				//
				j = 1;
				aux_a = ell_a[IDX(i, j)] * zor;
				aux_b = ell_b[IDX(i, j)];
				aux_c = ell_c[IDX(i, j)] * roz;
				aux_d = dz * ell_d[IDX(i, j)];
				aux_e = dr * ell_e[IDX(i, j)];
				aux_s = dr * dz * ell_s[IDX(i, j)];

				A.ia[IDX(i, j)] = BASE + offset;
				A.a[offset] = (-aux_a + aux_d) * twelfth;
				A.a[offset + 1] = (double)z_sym * (-aux_b / 48.0);
				A.a[offset + 2] = aux_b / 48.0;
				A.a[offset + 3] = aux_b * third;
				A.a[offset + 4] = 2.0 * (2.0 * aux_a - aux_d) * third;
				A.a[offset + 5] = -aux_b * third;
				A.a[offset + 6] = 2.0 * (2.0 * aux_c - aux_e) * third;
				A.a[offset + 7] = -2.5 * (aux_a + aux_c) + aux_s;
				A.a[offset + 8] = 2.0 * (2.0 * aux_c + aux_e) * third + (double)z_sym * ((-aux_c + aux_e) * twelfth);
				A.a[offset + 9] = -(aux_c + aux_e) * twelfth;
				A.a[offset + 10] = -aux_b * third;
				A.a[offset + 11] = 2.0 * (2.0 * aux_a + aux_d) * third;
				A.a[offset + 12] = aux_b * third;
				A.a[offset + 13] = -(aux_a + aux_d) * twelfth;
				A.a[offset + 14] = (double)z_sym * (aux_b / 48.0);
				A.a[offset + 15] = -aux_b / 48.0;
				A.ja[offset] = BASE + IDX(i - 2, j);
				A.ja[offset + 1] = BASE + IDX(i - 2, j + 1);
				A.ja[offset + 2] = BASE + IDX(i - 2, j + 2);
				A.ja[offset + 3] = BASE + IDX(i - 1, j - 1);
				A.ja[offset + 4] = BASE + IDX(i - 1, j);
				A.ja[offset + 5] = BASE + IDX(i - 1, j + 1);
				A.ja[offset + 6] = BASE + IDX(i, j - 1);
				A.ja[offset + 7] = BASE + IDX(i, j);
				A.ja[offset + 8] = BASE + IDX(i, j + 1);
				A.ja[offset + 9] = BASE + IDX(i, j + 2);
				A.ja[offset + 10] = BASE + IDX(i + 1, j - 1);
				A.ja[offset + 11] = BASE + IDX(i + 1, j);
				A.ja[offset + 12] = BASE + IDX(i + 1, j + 1);
				A.ja[offset + 13] = BASE + IDX(i + 2, j);
				A.ja[offset + 14] = BASE + IDX(i + 2, j + 1);
				A.ja[offset + 15] = BASE + IDX(i + 2, j + 2);
				ell_f[IDX(i, j)] *= dr * dz;
				offset += 16;

				// Centered interior points.
				//
				// o       o       o
				//   \     |     /
				//     o   o   o
				//       \ | /
				// o---o---x---o---o
				//       / | \
				//     o   o   o
				//   /     |     \
				// o       o       o
				//
				// Thus:
				// 
				// u(i-2, j-2) : -b/48
				// u(i-2, j  ) : (-a + d)/12
				// u(i-2, j+2) : b/48
				// u(i-1, j-1) : b/3
				// u(i-1, j  ) : 2(2a - d)/3
				// u(i-1, j+1) : -b/3
				// u(i  , j-2) : (-c + e)/12
				// u(i  , j-1) : 2(2c - e)/3
				// u(i  , j  ) : -5(a + c)/2 + s
				// u(i  , j+1) : 2(2c + e)/3
				// u(i  , j+2) : -(c + e)/12
				// u(i+1, j-1) : -b/3
				// u(i+1, j  ) : 2(2a + d)/3
				// u(i+1, j+1) : b/3
				// u(i+2, j-2) : b/48
				// u(i+2, j  ) : -(a + d)/12
				// u(i+2, j+2) : -b/48
				//
				for (j = 2; j < NzInterior; j++)
				{
					aux_a = ell_a[IDX(i, j)] * zor;
					aux_b = ell_b[IDX(i, j)];
					aux_c = ell_c[IDX(i, j)] * roz;
					aux_d = dz * ell_d[IDX(i, j)];
					aux_e = dr * ell_e[IDX(i, j)];
					aux_s = dr * dz * ell_s[IDX(i, j)];

					A.ia[IDX(i, j)] = BASE + offset;
					A.a[offset] = -aux_b / 48.0;
					A.a[offset + 1] = (-aux_a + aux_d) * twelfth;
					A.a[offset + 2] = aux_b / 48.0;
					A.a[offset + 3] = aux_b * third;
					A.a[offset + 4] = 2.0 * (2.0 * aux_a - aux_d) * third;
					A.a[offset + 5] = -aux_b * third;
					A.a[offset + 6] = (-aux_c + aux_e) * twelfth;
					A.a[offset + 7] = 2.0 * (2.0 * aux_c - aux_e) * third;
					A.a[offset + 8] = -2.5 * (aux_a + aux_c) + aux_s;
					A.a[offset + 9] = 2.0 * (2.0 * aux_c + aux_e) * third;
					A.a[offset + 10] = -(aux_c + aux_e) * twelfth;
					A.a[offset + 11] = -aux_b * third;
					A.a[offset + 12] = 2.0 * (2.0 * aux_a + aux_d) * third;
					A.a[offset + 13] = aux_b * third;
					A.a[offset + 14] = aux_b / 48.0;
					A.a[offset + 15] = -(aux_a + aux_d) * twelfth;
					A.a[offset + 16] = -aux_b / 48.0;
					A.ja[offset] = BASE + IDX(i - 2, j - 2);
					A.ja[offset + 1] = BASE + IDX(i - 2, j);
					A.ja[offset + 2] = BASE + IDX(i - 2, j + 2);
					A.ja[offset + 3] = BASE + IDX(i - 1, j - 1);
					A.ja[offset + 4] = BASE + IDX(i - 1, j);
					A.ja[offset + 5] = BASE + IDX(i - 1, j + 1);
					A.ja[offset + 6] = BASE + IDX(i, j - 2);
					A.ja[offset + 7] = BASE + IDX(i, j - 1);
					A.ja[offset + 8] = BASE + IDX(i, j);
					A.ja[offset + 9] = BASE + IDX(i, j + 1);
					A.ja[offset + 10] = BASE + IDX(i, j + 2);
					A.ja[offset + 11] = BASE + IDX(i + 1, j - 1);
					A.ja[offset + 12] = BASE + IDX(i + 1, j);
					A.ja[offset + 13] = BASE + IDX(i + 1, j + 1);
					A.ja[offset + 14] = BASE + IDX(i + 2, j - 2);
					A.ja[offset + 15] = BASE + IDX(i + 2, j);
					A.ja[offset + 16] = BASE + IDX(i + 2, j + 2);
					ell_f[IDX(i, j)] *= dr * dz;
					offset += 17;
				}

				// Stencil for top strip.
				//
				// o   o   o   o   o
				//   \   \ | /   /
				// o---o---x---o---o
				//   /   / | \   \
				// o   o   o   o   o
				//   /   / | \   \
				// o   o   o   o   o
				//   /   / | \   \
				// o   o   o   o   o
				//         |      
				//         o        
				//
				// Thus:
				// 
				// u(i-2, j-3) : -b/144
				// u(i-2, j-2) : b/24
				// u(i-2, j-1) : -b/8
				// u(i-2, j  ) : (-6a + 5b + 6d)/72
				// u(i-2, j+1) : b/48
				// u(i-1, j-3) : b/18
				// u(i-1, j-2) : -b/3
				// u(i-1, j-1) : b
				// u(i-1, j  ) : (12a - 5b - 6d)/9
				// u(i-1, j+1) : -b/6
				// u(i  , j-4) : c/12
				// u(i  , j-3) : -(6c + e)/12
				// u(i  , j-2) : (7c + 3e)/6
				// u(i  , j-1) : -(2c + 9e)/6
				// u(i  , j  ) : -5a/2 -5c/4 + 5e/6 + s
				// u(i  , j+1) : 5c/6 + e/4
				// u(i+1, j-3) : -b/18
				// u(i+1, j-2) : b/3
				// u(i+1, j-1) : -b
				// u(i+1, j  ) : (12a + 5b + 6d)/9
				// u(i+1, j+1) : b/6
				// u(i+2, j-3) : b/144
				// u(i+2, j-2) : -b/24
				// u(i+2, j-1) : b/8
				// u(i+2, j  ) : -(6a + 5b + 6d)/72
				// u(i+2, j+1) : -b/48
				//
				j = NzInterior;
				aux_a = ell_a[IDX(i, j)] * zor;
				aux_b = ell_b[IDX(i, j)];
				aux_c = ell_c[IDX(i, j)] * roz;
				aux_d = dz * ell_d[IDX(i, j)];
				aux_e = dr * ell_e[IDX(i, j)];
				aux_s = dr * dz * ell_s[IDX(i, j)];

				A.ia[IDX(i, j)] = BASE + offset;
				A.a[offset] = -aux_b / 144.0;
				A.a[offset + 1] = aux_b / 24.0;
				A.a[offset + 2] = -aux_b * 0.125;
				A.a[offset + 3] = (-6.0 * aux_a + 5.0 * aux_b + 6.0 * aux_d) / 72.0;
				A.a[offset + 4] = aux_b / 48.0;
				A.a[offset + 5] = aux_b / 18.0;
				A.a[offset + 6] = -aux_b * third;
				A.a[offset + 7] = aux_b;
				A.a[offset + 8] = (12.0 * aux_a - 5.0 * aux_b - 6.0 * aux_d) / 9.0;
				A.a[offset + 9] = -aux_b * sixth;
				A.a[offset + 10] = aux_c * twelfth;
				A.a[offset + 11] = -(6.0 * aux_c + aux_e) * twelfth;
				A.a[offset + 12] = (7.0 * aux_c + 3.0 * aux_e) * sixth;
				A.a[offset + 13] = -(2.0 * aux_c + 9.0 * aux_e) * sixth;
				A.a[offset + 14] = -2.5 * aux_a - 1.25 * aux_c + 5.0 * aux_e * sixth + aux_s;
				A.a[offset + 15] = 5.0 * aux_c * sixth + aux_e * 0.25;
				A.a[offset + 16] = -aux_b / 18.0;
				A.a[offset + 17] = aux_b * third;
				A.a[offset + 18] = -aux_b;
				A.a[offset + 19] = (12.0 * aux_a + 5.0 * aux_b + 6.0 * aux_d) / 9.0;
				A.a[offset + 20] = aux_b * sixth;
				A.a[offset + 21] = aux_b / 144.0;
				A.a[offset + 22] = -aux_b / 24.0;
				A.a[offset + 23] = aux_b * 0.125;
				A.a[offset + 24] = -(6.0 * aux_a + 5.0 * aux_b + 6.0 * aux_d) / 72.0;
				A.a[offset + 25] = -aux_b / 48.0;
				A.ja[offset] = BASE + IDX(i - 2, j - 3);
				A.ja[offset + 1] = BASE + IDX(i - 2, j - 2);
				A.ja[offset + 2] = BASE + IDX(i - 2, j - 1);
				A.ja[offset + 3] = BASE + IDX(i - 2, j);
				A.ja[offset + 4] = BASE + IDX(i - 2, j + 1);
				A.ja[offset + 5] = BASE + IDX(i - 1, j - 3);
				A.ja[offset + 6] = BASE + IDX(i - 1, j - 2);
				A.ja[offset + 7] = BASE + IDX(i - 1, j - 1);
				A.ja[offset + 8] = BASE + IDX(i - 1, j);
				A.ja[offset + 9] = BASE + IDX(i - 1, j + 1);
				A.ja[offset + 10] = BASE + IDX(i, j - 4);
				A.ja[offset + 11] = BASE + IDX(i, j - 3);
				A.ja[offset + 12] = BASE + IDX(i, j - 2);
				A.ja[offset + 13] = BASE + IDX(i, j - 1);
				A.ja[offset + 14] = BASE + IDX(i, j);
				A.ja[offset + 15] = BASE + IDX(i, j + 1);
				A.ja[offset + 16] = BASE + IDX(i + 1, j - 3);
				A.ja[offset + 17] = BASE + IDX(i + 1, j - 2);
				A.ja[offset + 18] = BASE + IDX(i + 1, j - 1);
				A.ja[offset + 19] = BASE + IDX(i + 1, j);
				A.ja[offset + 20] = BASE + IDX(i + 1, j + 1);
				A.ja[offset + 21] = BASE + IDX(i + 2, j - 3);
				A.ja[offset + 22] = BASE + IDX(i + 2, j - 2);
				A.ja[offset + 23] = BASE + IDX(i + 2, j - 1);
				A.ja[offset + 24] = BASE + IDX(i + 2, j);
				A.ja[offset + 25] = BASE + IDX(i + 2, j + 1);
				ell_f[IDX(i, j)] *= dr * dz;
				offset += 26;

				j = NzInterior + 1;
				z = (double)j - 0.5;
                		// Overall division by dz**2.
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
					case 2:
						robin2 = (rr2 / z) * (rr2 / z);
						robin1 = (rr2 / z) * (4.0 - (r * roz / z) * (r * roz / z));
						A.a[offset] = -5.0 * robin2 / 12.0;
						A.a[offset + 1] = 61.0 * robin2 / 24.0 + 0.125 * robin1;
						A.a[offset + 2] = -(6.5 * robin2 + 2.0 * robin1 / 3.0);
						A.a[offset + 3] = 107.0 * robin2 / 12.0 + 1.5 * robin1;
						A.a[offset + 4] = -(77.0 * robin2 / 12.0 + 2.0 * robin1);
						A.a[offset + 5] = 1.0 + 1.875 * robin2 + 25.0 * robin1 / 24.0;
						A.ja[offset] = BASE + IDX(i, NzInterior - 4);
						A.ja[offset + 1] = BASE + IDX(i, NzInterior - 3);
						A.ja[offset + 2] = BASE + IDX(i, NzInterior - 2);
						A.ja[offset + 3] = BASE + IDX(i, NzInterior - 1);
						A.ja[offset + 4] = BASE + IDX(i, NzInterior);
						A.ja[offset + 5] = BASE + IDX(i, NzInterior + 1);
						break;
					case 3:
						robin3 = (rr2 / z) * (rr2 / z) * (rr2 / z);
						robin2 = (rr2 / z) * (rr2 / z) * (9.0 - 3.0 * (r * roz / z) * (r * roz / z));
						robin1 = (rr2 / z) * (18.0 + (r * roz / z) * (r * roz / z) * (-9.0 + 3.0 * (rr2 / (z * z))));
						A.a[offset] = 0.3125 * robin3;
						A.a[offset + 1] = -(13.0 * robin3 / 6.0 + 5.0 * robin2 / 36.0);
						A.a[offset + 2] = 307.0 * robin3 / 48.0 + 61.0 * robin2 / 72.0 + robin1 / 24.0;
						A.a[offset + 3] = -(31.0 * robin3 / 3.0 + 13.0 * robin2 / 6.0 + 2.0 * robin1 / 9.0);
						A.a[offset + 4] = 461.0 * robin3 / 48.0 + 107.0 * robin2 / 36.0 + 0.5 * robin1;
						A.a[offset + 5] = -(29.0 * robin3 / 6.0 + 77.0 * robin2 / 36.0 + 2.0 * robin1 / 3.0);
						A.a[offset + 6] = 1.0 + 49.0 * robin3 / 48.0 + 0.625 * robin2 + 25.0 * robin1 / 72.0;
						A.ja[offset] = BASE + IDX(i, NzInterior - 5);
						A.ja[offset + 1] = BASE + IDX(i, NzInterior - 4);
						A.ja[offset + 2] = BASE + IDX(i, NzInterior - 3);
						A.ja[offset + 3] = BASE + IDX(i, NzInterior - 2);
						A.ja[offset + 4] = BASE + IDX(i, NzInterior - 1);
						A.ja[offset + 5] = BASE + IDX(i, NzInterior);
						A.ja[offset + 6] = BASE + IDX(i, NzInterior + 1);
						break;
				}
				ell_f[IDX(i, j)] = uInf;
			}
		}

#ifdef DEBUG
		printf("CSR FILL-IN: Stage 3.\n");
#endif

		// At this point we have filled:
		offset = (2 + 2) + 2 * NzInterior
			+ 2 + 14 + 16 * (NzInterior - 2) + 21 + n_robin
			+ (NrInterior - 2) * (2 + 16 + 17 * (NzInterior - 2) + 26 + n_robin);

		// offset = (2 + 2) + 2 * NzInterior 
		// + 2 + 14 + 16 * (NzInterior - 2) + 21 + n_robin.
		// + (NrInterior - 2) * (2 + 16 + 17 * (NzInterior - 2) + 26 + n_robin).

		// Interior bottom-right corner with equatorial symmetry.
		i = NrInterior;
		r = (double)i - 0.5;

		j = 0;
		A.ia[IDX(NrInterior, 0)] = BASE + offset;
		A.a[offset] = 1.0;
		A.a[offset + 1] = -(double)z_sym;
		A.ja[offset] = BASE + IDX(NrInterior, 0);
		A.ja[offset + 1] = BASE + IDX(NrInterior, 1);
		ell_f[IDX(NrInterior, 0)] = 0.0;
		offset += 2;

		// offset = (2 + 2) + 2 * NzInterior + 2 + 14 + 16 * (NzInterior - 2) + 21 + n_robin
		// 	+ (NrInterior - 2) * (2 + 16 + 17 * (NzInterior - 2) + 26 +  n_robin)
		// 	+ 2.


		// Right strip lower corner stencil.
		//
		//     o   o   o   o   o
		//       \   \   \ | /
		//     o   o   o   o   o
		//       \   \   \ | /
		// o---o---o---o---x---o
		//       /   /   / | \
		//     o   o   o   o   o
		// 
		// Thus:
		// 
		// u(i-4, j  ) : a/12
		// u(i-3, j-1) : b/18
		// u(i-3, j  ) : -(6a + d)/12
		// u(i-3, j+1) : -b/18 + z_sym*(-b/144)
		// u(i-3, j+2) : b/144
		// u(i-2, j-1) : -b/3
		// u(i-2, j  ) : (7a + 3d)/6
		// u(i-2, j+1) : b/3 + z_sym*(b/24)
		// u(i-2, j+2) : -b/24
		// u(i-1, j-1) : b
		// u(i-1, j  ) : -(2a + 9d)/6
		// u(i-1, j+1) : -b + z_sym*(-b/8)
		// u(i-1, j+2) : b/8
		// u(i  , j-1) : (-5b + 12c - 6e)/9
		// u(i  , j  ) : -5a/4 - 5c/2 + 5d/6 + s
		// u(i  , j+1) : (5b + 12c + 6e)/9 + z_sym*((5b - 6c + 6e)/72)
		// u(i  , j+2) : -(5b + 6c + 6e)/72
		// u(i+1, j-1) : -b/6
		// u(i+1, j  ) : 5a/6 + d/4
		// u(i+1, j+1) : b/6 + z_sym*(b/48)
		// u(i+1, j+2) : -b/48
		//
		j = 1;
		aux_a = ell_a[IDX(i, j)] * zor;
		aux_b = ell_b[IDX(i, j)];
		aux_c = ell_c[IDX(i, j)] * roz;
		aux_d = dz * ell_d[IDX(i, j)];
		aux_e = dr * ell_e[IDX(i, j)];
		aux_s = dr * dz * ell_s[IDX(i, j)];

		A.ia[IDX(i, j)] = BASE + offset;
		A.a[offset] = aux_a * twelfth;
		A.a[offset + 1] = aux_b / 18.0;
		A.a[offset + 2] = -(6.0 * aux_a + aux_d) * twelfth;
		A.a[offset + 3] = -aux_b / 18.0 + (double)z_sym * (-aux_b / 144.0);
		A.a[offset + 4] = aux_b / 144.0;
		A.a[offset + 5] = -aux_b * third;
		A.a[offset + 6] = (7.0 * aux_a + 3.0 * aux_d) * sixth;
		A.a[offset + 7] = aux_b * third + (double)z_sym * (aux_b / 24.0);
		A.a[offset + 8] = -aux_b / 24.0;
		A.a[offset + 9] = aux_b;
		A.a[offset + 10] = -(2.0 * aux_a + 9.0 * aux_d) * sixth;
		A.a[offset + 11] = -aux_b + (double)z_sym * (-aux_b * 0.125);
		A.a[offset + 12] = aux_b * 0.125;
		A.a[offset + 13] = (-5.0 * aux_b + 12.0 * aux_c - 6.0 * aux_e) / 9.0;
		A.a[offset + 14] = -1.25 * aux_a - 2.5 * aux_c + 5.0 * aux_d * sixth + aux_s;
		A.a[offset + 15] = (5.0 * aux_b + 12.0 * aux_c + 6.0 * aux_e) / 9.0 + (double)z_sym * ((5.0 * aux_b - 6.0 * aux_c + 6.0 * aux_e) / 72.0);
		A.a[offset + 16] = -(5.0 * aux_b + 6.0 * aux_c + 6.0 * aux_e) / 72.0;
		A.a[offset + 17] = -aux_b * sixth;
		A.a[offset + 18] = 5.0 * aux_a * sixth + aux_d * 0.25;
		A.a[offset + 19] = aux_b * sixth + (double)z_sym * (aux_b / 48.0);
		A.a[offset + 20] = -aux_b / 48.0;
		A.ja[offset] = BASE + IDX(i - 4, j);
		A.ja[offset + 1] = BASE + IDX(i - 3, j - 1);
		A.ja[offset + 2] = BASE + IDX(i - 3, j);
		A.ja[offset + 3] = BASE + IDX(i - 3, j + 1);
		A.ja[offset + 4] = BASE + IDX(i - 3, j + 2);
		A.ja[offset + 5] = BASE + IDX(i - 2, j - 1);
		A.ja[offset + 6] = BASE + IDX(i - 2, j);
		A.ja[offset + 7] = BASE + IDX(i - 2, j + 1);
		A.ja[offset + 8] = BASE + IDX(i - 2, j + 2);
		A.ja[offset + 9] = BASE + IDX(i - 1, j - 1);
		A.ja[offset + 10] = BASE + IDX(i - 1, j);
		A.ja[offset + 11] = BASE + IDX(i - 1, j + 1);
		A.ja[offset + 12] = BASE + IDX(i - 1, j + 2);
		A.ja[offset + 13] = BASE + IDX(i, j - 1);
		A.ja[offset + 14] = BASE + IDX(i, j);
		A.ja[offset + 15] = BASE + IDX(i, j + 1);
		A.ja[offset + 16] = BASE + IDX(i, j + 2);
		A.ja[offset + 17] = BASE + IDX(i + 1, j - 1);
		A.ja[offset + 18] = BASE + IDX(i + 1, j);
		A.ja[offset + 19] = BASE + IDX(i + 1, j + 1);
		A.ja[offset + 20] = BASE + IDX(i + 1, j + 2);
		ell_f[IDX(i, j)] *= dr * dz;
		offset += 21;

		// offset = (2 + 2) + 2 * NzInterior + 2 + 14 + 16 * (NzInterior - 2) + 21 + robin
		// 	+ (NrInterior - 2) * (2 + 16 + 17 * (NzInterior - 2) + 26 +  n_robin)
		// 	+ 2 + 21.

		// Set temporary offset.
		t_offset = offset;

        	#pragma omp parallel shared(A, ell_f) private(offset,\
		aux_a, aux_b, aux_c, aux_d, aux_e, aux_s)
		{
            		#pragma omp for schedule(guided)
			for (j = 2; j < NzInterior; j++)
			{
				// Eac iteration fills 26 elements.
				offset = t_offset + (j - 2) * 26;

				// Right strip stencil.
				//
				//     o   o   o   o   o
				//       \   \   \ | /
				//     o   o   o   o   o
				//       \   \   \ | /
				// o---o---o---o---x---o
				//       /   /   / | \
				//     o   o   o   o   o
				//       /   /   / | \
				//     o   o   o   o   o
				// 
				// Thus:
				// 
				// u(i-4, j  ) : a/12
				// u(i-3, j-2) : -b/144
				// u(i-3, j-1) : b/18
				// u(i-3, j  ) : -(6a + d)/12
				// u(i-3, j+1) : -b/18
				// u(i-3, j+2) : b/144
				// u(i-2, j-2) : b/24
				// u(i-2, j-1) : -b/3
				// u(i-2, j  ) : (7a + 3d)/6
				// u(i-2, j+1) : b/3
				// u(i-2, j+2) : -b/24
				// u(i-1, j-2) : -b/8
				// u(i-1, j-1) : b
				// u(i-1, j  ) : -(2a + 9d)/6
				// u(i-1, j+1) : -b
				// u(i-1, j+2) : b/8
				// u(i  , j-2) : (5b - 6c + 6e)/72
				// u(i  , j-1) : (-5b + 12c - 6e)/9
				// u(i  , j  ) : -5a/4 - 5c/2 + 5d/6 + s
				// u(i  , j+1) : (5b + 12c + 6e)/9
				// u(i  , j+2) : -(5b + 6c + 6e)/72
				// u(i+1, j-2) : b/48
				// u(i+1, j-1) : -b/6
				// u(i+1, j  ) : 5a/6 + d/4
				// u(i+1, j+1) : b/6
				// u(i+1, j+2) : -b/48
				//
				aux_a = ell_a[IDX(i, j)] * zor;
				aux_b = ell_b[IDX(i, j)];
				aux_c = ell_c[IDX(i, j)] * roz;
				aux_d = dz * ell_d[IDX(i, j)];
				aux_e = dr * ell_e[IDX(i, j)];
				aux_s = dr * dz * ell_s[IDX(i, j)];

				A.ia[IDX(i, j)] = BASE + offset;
				A.a[offset] = aux_a * twelfth;
				A.a[offset + 1] = -aux_b / 144.0;
				A.a[offset + 2] = aux_b / 18.0;
				A.a[offset + 3] = -(6.0 * aux_a + aux_d) * twelfth;
				A.a[offset + 4] = -aux_b / 18.0;
				A.a[offset + 5] = aux_b / 144.0;
				A.a[offset + 6] = aux_b / 24.0;
				A.a[offset + 7] = -aux_b * third;
				A.a[offset + 8] = (7.0 * aux_a + 3.0 * aux_d) * sixth;
				A.a[offset + 9] = aux_b * third;
				A.a[offset + 10] = -aux_b / 24.0;
				A.a[offset + 11] = -aux_b * 0.125;
				A.a[offset + 12] = aux_b;
				A.a[offset + 13] = -(2.0 * aux_a + 9.0 * aux_d) * sixth;
				A.a[offset + 14] = -aux_b;
				A.a[offset + 15] = aux_b * 0.125;
				A.a[offset + 16] = (5.0 * aux_b - 6.0 * aux_c + 6.0 * aux_e) / 72.0;
				A.a[offset + 17] = (-5.0 * aux_b + 12.0 * aux_c - 6.0 * aux_e) / 9.0;
				A.a[offset + 18] = -1.25 * aux_a - 2.5 * aux_c + 5.0 * aux_d * sixth + aux_s;
				A.a[offset + 19] = (5.0 * aux_b + 12.0 * aux_c + 6.0 * aux_e) / 9.0;
				A.a[offset + 20] = -(5.0 * aux_b + 6.0 * aux_c + 6.0 * aux_e) / 72.0;
				A.a[offset + 21] = aux_b / 48.0;
				A.a[offset + 22] = -aux_b * sixth;
				A.a[offset + 23] = 5.0 * aux_a * sixth + aux_d * 0.25;
				A.a[offset + 24] = aux_b * sixth;
				A.a[offset + 25] = -aux_b / 48.0;
				A.ja[offset] = BASE + IDX(i - 4, j);
				A.ja[offset + 1] = BASE + IDX(i - 3, j - 2);
				A.ja[offset + 2] = BASE + IDX(i - 3, j - 1);
				A.ja[offset + 3] = BASE + IDX(i - 3, j);
				A.ja[offset + 4] = BASE + IDX(i - 3, j + 1);
				A.ja[offset + 5] = BASE + IDX(i - 3, j + 2);
				A.ja[offset + 6] = BASE + IDX(i - 2, j - 2);
				A.ja[offset + 7] = BASE + IDX(i - 2, j - 1);
				A.ja[offset + 8] = BASE + IDX(i - 2, j);
				A.ja[offset + 9] = BASE + IDX(i - 2, j + 1);
				A.ja[offset + 10] = BASE + IDX(i - 2, j + 2);
				A.ja[offset + 11] = BASE + IDX(i - 1, j - 2);
				A.ja[offset + 12] = BASE + IDX(i - 1, j - 1);
				A.ja[offset + 13] = BASE + IDX(i - 1, j);
				A.ja[offset + 14] = BASE + IDX(i - 1, j + 1);
				A.ja[offset + 15] = BASE + IDX(i - 1, j + 2);
				A.ja[offset + 16] = BASE + IDX(i, j - 2);
				A.ja[offset + 17] = BASE + IDX(i, j - 1);
				A.ja[offset + 18] = BASE + IDX(i, j);
				A.ja[offset + 19] = BASE + IDX(i, j + 1);
				A.ja[offset + 20] = BASE + IDX(i, j + 2);
				A.ja[offset + 21] = BASE + IDX(i + 1, j - 2);
				A.ja[offset + 22] = BASE + IDX(i + 1, j - 1);
				A.ja[offset + 23] = BASE + IDX(i + 1, j);
				A.ja[offset + 24] = BASE + IDX(i + 1, j + 1);
				A.ja[offset + 25] = BASE + IDX(i + 1, j + 2);
				ell_f[IDX(i, j)] *= dr * dz;
			}
		}


#ifdef DEBUG
		printf("CSR FILL-IN: Stage 4.\n");
#endif

		// At this point we have filled:
		offset = (2 + 2) + 2 * NzInterior + 2 + 14 + 16 * (NzInterior - 2) + 21 + n_robin
			+ (NrInterior - 2) * (2 + 16 + 17 * (NzInterior - 2) + 26 + n_robin)
			+ 2 + 21 + 26 * (NzInterior - 2);

		// offset = (2 + 2) + 2 * NzInterior + 2 + 14 + 16 * (NzInterior - 2) + 21 + n_robin
		// 	+ (NrInterior - 2) * (2 + 16 + 17 * (NzInterior - 2) + 26 +  n_robin)
		// 	+ 2 + 21 + 26 * (NzInterior - 2).

		// Fill interior top-right corner.
		//
		//     o   o   o   o   o
		//       \   \   \ | /
		// o---o---o---o---x---o
		//       \  \    / | \
		//     o   o   o   o   o
		//       \   /   \ | \
		//     o   o   o   o   o
		//       /   \   \ | \
		//     o   o   o   o   o
		//                 |
		//                 o   
		//
		// Thus:
		//
		// u(i-4, j  ) : a/12
		// u(i-3, j-3) : b/144
		// u(i-3, j-2) : -b/24
		// u(i-3, j-1) : b/8
		// u(i-3, j  ) : -(36a + 5b + 6d)/72
		// u(i-3, j+1) : -b/48
		// u(i-2, j-3) : -b/24
		// u(i-2, j-2) : b/4
		// u(i-2, j-1) : -3b/4
		// u(i-2, j  ) : (14a + 5b + 6d)/12
		// u(i-2, j+1) : b/8
		// u(i-1, j-3) : b/8
		// u(i-1, j-2) : -3b/4
		// u(i-1, j-1) : 9b/4
		// u(i-1, j  ) : -a/3 - 5b/4 - 3d/2
		// u(i-1, j+1) : -3b/8
		// u(i  , j-4) : c/12
		// u(i  , j-3) : -(36c + 5b + 6e)/72
		// u(i  , j-2) : (14c + 5b + 6e)/12
		// u(i  , j-1) : -c/3 - 5b/4 - 3e/2
		// u(i  , j  ) : -5(9a - 5b + 9c - 6d - 6e)/36 + s
		// u(i  , j+1) : (5b + 20c + 6e)/24
		// u(i+1, j-3) : -b/48
		// u(i+1, j-2) : b/8
		// u(i+1, j-1) : -3b/8
		// u(i+1, j  ) : (20a + 5b + 6d)/24
		// u(i+1, j+1) : b/16
		//
		j = NzInterior;
		aux_a = ell_a[IDX(i, j)] * zor;
		aux_b = ell_b[IDX(i, j)];
		aux_c = ell_c[IDX(i, j)] * roz;
		aux_d = dz * ell_d[IDX(i, j)];
		aux_e = dr * ell_e[IDX(i, j)];
		aux_s = dr * dz * ell_s[IDX(i, j)];

		A.ia[IDX(i, j)] = BASE + offset;
		A.a[offset] = aux_a * twelfth;
		A.a[offset + 1] = aux_b / 144.0;
		A.a[offset + 2] = -aux_b / 24.0;
		A.a[offset + 3] = aux_b * 0.125;
		A.a[offset + 4] = -(36.0 * aux_a + 5.0 * aux_b + 6.0 * aux_d) / 72.0;
		A.a[offset + 5] = -aux_b / 48.0;
		A.a[offset + 6] = -aux_b / 24.0;
		A.a[offset + 7] = aux_b * 0.25;
		A.a[offset + 8] = -aux_b * 0.75;
		A.a[offset + 9] = (14.0 * aux_a + 5.0 * aux_b + 6.0 * aux_d) * twelfth;
		A.a[offset + 10] = aux_b * 0.125;
		A.a[offset + 11] = aux_b * 0.125;
		A.a[offset + 12] = -aux_b * 0.75;
		A.a[offset + 13] = 9.0 * aux_b * 0.25;
		A.a[offset + 14] = -aux_a * third - 5.0 * aux_b * 0.25 - 1.5 * aux_d;
		A.a[offset + 15] = -3.0 * aux_b * 0.125;
		A.a[offset + 16] = aux_c * twelfth;
		A.a[offset + 17] = -(36.0 * aux_c + 5.0 * aux_b + 6.0 * aux_e) / 72.0;
		A.a[offset + 18] = (14.0 * aux_c + 5.0 * aux_b + 6.0 * aux_e) * twelfth;
		A.a[offset + 19] = -aux_c * third - 5.0 * aux_b * 0.25 - 1.5 * aux_e;
		A.a[offset + 20] = -5.0 * (9.0 * aux_a - 5.0 * aux_b + 9.0 * aux_c - 6.0 * aux_d - 6.0 * aux_e) / 36.0 + aux_s;
		A.a[offset + 21] = (5.0 * aux_b + 20.0 * aux_c + 6.0 * aux_e) / 24.0;
		A.a[offset + 22] = -aux_b / 48.0;
		A.a[offset + 23] = aux_b * 0.125;
		A.a[offset + 24] = -3.0 * aux_b * 0.125;
		A.a[offset + 25] = (5.0 * aux_b + 20.0 * aux_a + 6.0 * aux_d) / 24.0;
		A.a[offset + 26] = aux_b * 0.0625;
		A.ja[offset] = BASE + IDX(i - 4, j);
		A.ja[offset + 1] = BASE + IDX(i - 3, j - 3);
		A.ja[offset + 2] = BASE + IDX(i - 3, j - 2);
		A.ja[offset + 3] = BASE + IDX(i - 3, j - 1);
		A.ja[offset + 4] = BASE + IDX(i - 3, j);
		A.ja[offset + 5] = BASE + IDX(i - 3, j + 1);
		A.ja[offset + 6] = BASE + IDX(i - 2, j - 3);
		A.ja[offset + 7] = BASE + IDX(i - 2, j - 2);
		A.ja[offset + 8] = BASE + IDX(i - 2, j - 1);
		A.ja[offset + 9] = BASE + IDX(i - 2, j);
		A.ja[offset + 10] = BASE + IDX(i - 2, j + 1);
		A.ja[offset + 11] = BASE + IDX(i - 1, j - 3);
		A.ja[offset + 12] = BASE + IDX(i - 1, j - 2);
		A.ja[offset + 13] = BASE + IDX(i - 1, j - 1);
		A.ja[offset + 14] = BASE + IDX(i - 1, j);
		A.ja[offset + 15] = BASE + IDX(i - 1, j + 1);
		A.ja[offset + 16] = BASE + IDX(i, j - 4);
		A.ja[offset + 17] = BASE + IDX(i, j - 3);
		A.ja[offset + 18] = BASE + IDX(i, j - 2);
		A.ja[offset + 19] = BASE + IDX(i, j - 1);
		A.ja[offset + 20] = BASE + IDX(i, j);
		A.ja[offset + 21] = BASE + IDX(i, j + 1);
		A.ja[offset + 22] = BASE + IDX(i + 1, j - 3);
		A.ja[offset + 23] = BASE + IDX(i + 1, j - 2);
		A.ja[offset + 24] = BASE + IDX(i + 1, j - 1);
		A.ja[offset + 25] = BASE + IDX(i + 1, j);
		A.ja[offset + 26] = BASE + IDX(i + 1, j + 1);
		ell_f[IDX(i, j)] *= dr * dz;
		offset += 27;

		// offset = (2 + 2) + 2 * NzInterior + 2 + 14 + 16 * (NzInterior - 2) + 21 + n_robin
		// 	+ (NrInterior - 2) * (2 + 16 + 17 * (NzInterior - 2) + 26 +  n_robin)
		// 	+ 2 + 21 + 26 * (NzInterior - 2) + 27

		j = NzInterior + 1;
		z = (double)j - 0.5;
        	// Overall division by dz**2.
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
			case 2:
				robin2 = (rr2 / z) * (rr2 / z);
				robin1 = (rr2 / z) * (4.0 - (r * roz / z) * (r * roz / z));
				A.a[offset] = -5.0 * robin2 / 12.0;
				A.a[offset + 1] = 61.0 * robin2 / 24.0 + 0.125 * robin1;
				A.a[offset + 2] = -(6.5 * robin2 + 2.0 * robin1 / 3.0);
				A.a[offset + 3] = 107.0 * robin2 / 12.0 + 1.5 * robin1;
				A.a[offset + 4] = -(77.0 * robin2 / 12.0 + 2.0 * robin1);
				A.a[offset + 5] = 1.0 + 1.875 * robin2 + 25.0 * robin1 / 24.0;
				A.ja[offset] = BASE + IDX(i, NzInterior - 4);
				A.ja[offset + 1] = BASE + IDX(i, NzInterior - 3);
				A.ja[offset + 2] = BASE + IDX(i, NzInterior - 2);
				A.ja[offset + 3] = BASE + IDX(i, NzInterior - 1);
				A.ja[offset + 4] = BASE + IDX(i, NzInterior);
				A.ja[offset + 5] = BASE + IDX(i, NzInterior + 1);
				break;
			case 3:
				robin3 = (rr2 / z) * (rr2 / z) * (rr2 / z);
				robin2 = (rr2 / z) * (rr2 / z) * (9.0 - 3.0 * (r * roz / z) * (r * roz / z));
				robin1 = (rr2 / z) * (18.0 + (r * roz / z) * (r * roz / z) * (-9.0 + 3.0 * (rr2 / (z * z))));
				A.a[offset] = 0.3125 * robin3;
				A.a[offset + 1] = -(13.0 * robin3 / 6.0 + 5.0 * robin2 / 36.0);
				A.a[offset + 2] = 307.0 * robin3 / 48.0 + 61.0 * robin2 / 72.0 + robin1 / 24.0;
				A.a[offset + 3] = -(31.0 * robin3 / 3.0 + 13.0 * robin2 / 6.0 + 2.0 * robin1 / 9.0);
				A.a[offset + 4] = 461.0 * robin3 / 48.0 + 107.0 * robin2 / 36.0 + 0.5 * robin1;
				A.a[offset + 5] = -(29.0 * robin3 / 6.0 + 77.0 * robin2 / 36.0 + 2.0 * robin1 / 3.0);
				A.a[offset + 6] = 1.0 + 49.0 * robin3 / 48.0 + 0.625 * robin2 + 25.0 * robin1 / 72.0;
				A.ja[offset] = BASE + IDX(i, NzInterior - 5);
				A.ja[offset + 1] = BASE + IDX(i, NzInterior - 4);
				A.ja[offset + 2] = BASE + IDX(i, NzInterior - 3);
				A.ja[offset + 3] = BASE + IDX(i, NzInterior - 2);
				A.ja[offset + 4] = BASE + IDX(i, NzInterior - 1);
				A.ja[offset + 5] = BASE + IDX(i, NzInterior);
				A.ja[offset + 6] = BASE + IDX(i, NzInterior + 1);
				break;

		}
		ell_f[IDX(i, j)] = uInf;
		offset += n_robin;

		// offset = (2 + 2) + 2 * NzInterior + 2 + 14 + 16 * (NzInterior - 2) + 21 + n_robin
		// 	+ (NrInterior - 2) * (2 + 16 + 17 * (NzInterior - 2) + 26 +  n_robin)
		// 	+ 2 + 21 + 26 * (NzInterior - 2) + 27 + n_robin.

		// Do bottom-right corner.
		i = NrInterior + 1;
		r = (double)i - 0.5;

		j = 0;
		A.ia[IDX(NrInterior + 1, 0)] = BASE + offset;
		A.a[offset] = 1.0;
		A.a[offset + 1] = -(double)z_sym;
		A.ja[offset] = BASE + IDX(NrInterior + 1, 0);
		A.ja[offset + 1] = BASE + IDX(NrInterior + 1, 1);
		ell_f[IDX(NrInterior + 1, 0)] = 0.0;
		offset += 2;

		// offset = (2 + 2 + 2) + 2 * NzInterior + 2 + 14 + 16 * (NzInterior - 2) + 21 + n_robin
		// 	+ (NrInterior - 2) * (2 + 16 + 17 * (NzInterior - 2) + 26 +  n_robin)
		// 	+ 2 + 21 + 26 * (NzInterior - 2) + 27 + n_robin.

		// Set temporary offset.
		t_offset = offset;

        	#pragma omp parallel shared(A, ell_f) private(offset, z, rr2,\
		robin1, robin2, robin3)
		{
			#pragma omp for schedule(guided)
			for (j = 1; j < NzInterior + 1; j++)
			{
				// Each iteration fills n_robin elements.
				offset = t_offset + n_robin * (j - 1);

				// Z coordinate.
				z = (double)j - 0.5;
                		// Overall division by dr**2.
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
					case 2:
						robin2 = (rr2 / r) * (rr2 / r);
						robin1 = (rr2 / r) * (4.0 - (z * zor / r) * (z * zor / r));
						A.a[offset] = -5.0 * robin2 / 12.0;
						A.a[offset + 1] = 61.0 * robin2 / 24.0 + 0.125 * robin1;
						A.a[offset + 2] = -(6.5 * robin2 + 2.0 * robin1 / 3.0);
						A.a[offset + 3] = 107.0 * robin2 / 12.0 + 1.5 * robin1;
						A.a[offset + 4] = -(77.0 * robin2 / 12.0 + 2.0 * robin1);
						A.a[offset + 5] = 1.0 + 1.875 * robin2 + 25.0 * robin1 / 24.0;
						A.ja[offset] = BASE + IDX(NrInterior - 4, j);
						A.ja[offset + 1] = BASE + IDX(NrInterior - 3, j);
						A.ja[offset + 2] = BASE + IDX(NrInterior - 2, j);
						A.ja[offset + 3] = BASE + IDX(NrInterior - 1, j);
						A.ja[offset + 4] = BASE + IDX(NrInterior, j);
						A.ja[offset + 5] = BASE + IDX(NrInterior + 1, j);
						break;
					case 3:
						robin3 = (rr2 / r) * (rr2 / r) * (rr2 / r);
						robin2 = (rr2 / r) * (rr2 / r) * (9.0 - 3.0 * (z * zor / r) * (z * zor / r));
						robin1 = (rr2 / r) * (18.0 + (z * zor / r) * (z * zor / r) * (-9.0 + 3.0 * (rr2 / (r * r))));
						A.a[offset] = 0.3125 * robin3;
						A.a[offset + 1] = -(13.0 * robin3 / 6.0 + 5.0 * robin2 / 36.0);
						A.a[offset + 2] = 307.0 * robin3 / 48.0 + 61.0 * robin2 / 72.0 + robin1 / 24.0;
						A.a[offset + 3] = -(31.0 * robin3 / 3.0 + 13.0 * robin2 / 6.0 + 2.0 * robin1 / 9.0);
						A.a[offset + 4] = 461.0 * robin3 / 48.0 + 107.0 * robin2 / 36.0 + 0.5 * robin1;
						A.a[offset + 5] = -(29.0 * robin3 / 6.0 + 77.0 * robin2 / 36.0 + 2.0 * robin1 / 3.0);
						A.a[offset + 6] = 1.0 + 49.0 * robin3 / 48.0 + 0.625 * robin2 + 25.0 * robin1 / 72.0;
						A.ja[offset] = BASE + IDX(NrInterior - 5, j);
						A.ja[offset + 1] = BASE + IDX(NrInterior - 4, j);
						A.ja[offset + 2] = BASE + IDX(NrInterior - 3, j);
						A.ja[offset + 3] = BASE + IDX(NrInterior - 2, j);
						A.ja[offset + 4] = BASE + IDX(NrInterior - 1, j);
						A.ja[offset + 5] = BASE + IDX(NrInterior, j);
						A.ja[offset + 6] = BASE + IDX(NrInterior + 1, j);
						break;
				}
				ell_f[IDX(i, j)] = uInf;
			}
		}

#ifdef DEBUG
		printf("CSR FILL-IN: Stage 5.\n");
#endif

		// At this point we have filled:
		offset = (2 + 2 + 2) + 2 * NzInterior + 2 + 14 + 16 * (NzInterior - 2) + 21 + n_robin
			+ (NrInterior - 2) * (2 + 16 + 17 * (NzInterior - 2) + 26 + n_robin)
			+ 2 + 21 + 26 * (NzInterior - 2) + 27 + n_robin
			+ n_robin * NzInterior;

		// offset = (2 + 2 + 2) + 2 * NzInterior + 2 + 14 + 16 * (NzInterior - 2) + 21 + n_robin
		// 	+ (NrInterior - 2) * (2 + 16 + 17 * (NzInterior - 2) + 26 +  n_robin)
		// 	+ 2 + 21 + 26 * (NzInterior - 2) + 27 + n_robin.
		// 	+ n_robin * NzInterior.

		// Finally, fill top corner with Robin
		j = NzInterior + 1;
		z = (double)j - 0.5;
        	// Radial coordinate over radial step.
		rrodrr = sqrt((r * r * dr * dr + z * z * dz * dz) / (dr * dr + dz * dz));
		A.ia[IDX(i, j)] = BASE + offset;
		switch (robin)
		{
			case 1:
				robin1 = rrodrr;
				A.a[offset] = robin1 / 4.0;
				A.a[offset + 1] = -4.0 * robin1 / 3.0;
				A.a[offset + 2] = 3.0 * robin1;
				A.a[offset + 3] = -4.0 * robin1;
				A.a[offset + 4] = 1.0 + 25.0 * robin1 / 12.0;
				A.ja[offset] = BASE + IDX(NrInterior - 3, NzInterior - 3);
				A.ja[offset + 1] = BASE + IDX(NrInterior - 2, NzInterior - 2);
				A.ja[offset + 2] = BASE + IDX(NrInterior - 1, NzInterior - 1);
				A.ja[offset + 3] = BASE + IDX(NrInterior, NzInterior);
				A.ja[offset + 4] = BASE + IDX(NrInterior + 1, NzInterior + 1);
				break;
			case 2:
				robin2 = rrodrr * rrodrr;
				robin1 = 4.0 * rrodrr;
				A.a[offset] = -5.0 * robin2 / 12.0;
				A.a[offset + 1] = 61.0 * robin2 / 24.0 + 0.125 * robin1;
				A.a[offset + 2] = -(6.5 * robin2 + 2.0 * robin1 / 3.0);
				A.a[offset + 3] = 107.0 * robin2 / 12.0 + 1.5 * robin1;
				A.a[offset + 4] = -(77.0 * robin2 / 12.0 + 2.0 * robin1);
				A.a[offset + 5] = 1.0 + 1.875 * robin2 + 25.0 * robin1 / 24.0;
				A.ja[offset] = BASE + IDX(NrInterior - 4, NzInterior - 4);
				A.ja[offset + 1] = BASE + IDX(NrInterior - 3, NzInterior - 3);
				A.ja[offset + 2] = BASE + IDX(NrInterior - 2, NzInterior - 2);
				A.ja[offset + 3] = BASE + IDX(NrInterior - 1, NzInterior - 1);
				A.ja[offset + 4] = BASE + IDX(NrInterior, NzInterior);
				A.ja[offset + 5] = BASE + IDX(NrInterior + 1, NzInterior + 1);
				break;
			case 3:
				robin3 = rrodrr * rrodrr * rrodrr;
				robin2 = 9.0 * rrodrr * rrodrr;
				robin1 = 18.0 * rrodrr;
				A.a[offset] = 0.3125 * robin3;
				A.a[offset + 1] = -(13.0 * robin3 / 6.0 + 5.0 * robin2 / 36.0);
				A.a[offset + 2] = 307.0 * robin3 / 48.0 + 61.0 * robin2 / 72.0 + robin1 / 24.0;
				A.a[offset + 3] = -(31.0 * robin3 / 3.0 + 13.0 * robin2 / 6.0 + 2.0 * robin1 / 9.0);
				A.a[offset + 4] = 461.0 * robin3 / 48.0 + 107.0 * robin2 / 36.0 + 0.5 * robin1;
				A.a[offset + 5] = -(29.0 * robin3 / 6.0 + 77.0 * robin2 / 36.0 + 2.0 * robin1 / 3.0);
				A.a[offset + 6] = 1.0 + 49.0 * robin3 / 48.0 + 0.625 * robin2 + 25.0 * robin1 / 72.0;
				A.ja[offset] = BASE + IDX(NrInterior - 5, NzInterior - 5);
				A.ja[offset + 1] = BASE + IDX(NrInterior - 4, NzInterior - 4);
				A.ja[offset + 2] = BASE + IDX(NrInterior - 3, NzInterior - 3);
				A.ja[offset + 3] = BASE + IDX(NrInterior - 2, NzInterior - 2);
				A.ja[offset + 4] = BASE + IDX(NrInterior - 1, NzInterior - 1);
				A.ja[offset + 5] = BASE + IDX(NrInterior, NzInterior);
				A.ja[offset + 6] = BASE + IDX(NrInterior + 1, NzInterior + 1);
				break;

		}
		ell_f[IDX(i, j)] = uInf;
		offset += n_robin;

		// offset = (2 + 2 + 2 + n_robin) + 2 * NzInterior 
		// + 2 + 14 + 16 * (NzInterior - 2) + 21 + n_robin
		// + (NrInterior - 2) * (2 + 16 + 17 * (NzInterior - 2) + 26 + n_robin)
		// + 2 + 21 + 26 * (NzInterior - 2) + 27 + n_robin.
		// + n_robin * NzInterior.

		// Assert fill-in.
		assert(offset == nnz);

		// Fill last element.
		A.ia[IDX(NrInterior + 1, NzInterior + 1) + 1] = BASE + nnz;
	}

#ifdef DEBUG
	csr_print(A, "ge_A_a.asc", "ge_A_ia.asc", "ge_A_ja.asc");
#endif
	// All done.
	return;
}
