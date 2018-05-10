// Global header for tools.
#include "tools.h"

// Reduce array u to elliptic solver-sized array g_u.
// 
// The number of ghost zones is not changed during execution.
void ghost_reduce(const double *u, 		// Original array to reduce.
			double *g_u,		// Output reduced array. 
			const int NrInterior, 	// Number of interior points in r.
			const int NzInterior, 	// Number of interior points in z.
			const int ghost)	// Number of ghost zones.
{
	// Auxiliary integers.
	int i, j;
	int NrTotal = ghost + NrInterior + 1;
	int NzTotal = ghost + NzInterior + 1;
	// Set this integer to the ghost zone just left(below) the r(z) axis.
	int k = ghost - 1;

	// Loop over interior points.
	#pragma omp parallel shared(g_u) private(j)
	{
		#pragma omp for schedule(guided)
		for (i = 0; i < NrInterior + 2; i++)
		{
			for (j = 0; j < NzInterior + 2; j++)
			{
				// g_u has a single ghost zone.
				g_u[i * (NzInterior + 2) + j] = u[IDX(k + i, k + j)];
			}
		}
	}

	return;
}

// Fill and array u using an elliptic solver-sized array g_u using symmetry conditions.
// 
// The number of ghost zones is not changed during execution. Therefore, ghost has to
// be reset to its original value.
void ghost_fill(const double *g_u, 	// Reduced to array to extend.
		double *u, 		// Output extended array.
		const int r_sym, 	// R symmetry condition: 1(even), -1(odd).
		const int z_sym, 	// Z symmetry condition: 1(even), -1(odd).
		const int NrInterior, 	// Number of interior points in r.
		const int NzInterior, 	// Number of interior points in z.
		const int ghost)	// Number of ghost zones.
{
	// Auxiliary integers.
	int i, j;
	int NrTotal = ghost + NrInterior + 1;
	int NzTotal = ghost + NzInterior + 1;
	// Set this integer to the ghost zone just left(below) the r(z) axis.
	int k = ghost - 1;

	// Fill points that coincide with reduced array.
	#pragma omp parallel shared(u) private(j)
	{
		#pragma omp for schedule(guided)
		for (i = 0; i < NrInterior + 2; i++)
		{
			for (j = 0; j < NzInterior + 2; j++)
			{
				// g_u has a single ghost zone that coincindes with ghost zone k.
				u[IDX(k + i, k + j)] = g_u[i * (NzInterior + 2) + j];
			}
		}
	}

	// Now fill remaining ghost zones:
	for (k = ghost - 2; k >= 0; k--)
	{
		// Correct R boundaries.
		#pragma omp parallel shared(u) 
		{
			#pragma omp for schedule(guided)
			for (j = k + 1; j < NzInterior + ghost + 1; j++)
			{
				// Symmetry.
				u[IDX(k, j)] = (double)r_sym * u[IDX(2 * ghost - 1 - k, j)];
			}
		}

		// Correct Z boundaries.
		#pragma omp parallel shared(u)
		{
			#pragma omp for schedule(guided)
			for (i = k + 1; i < NrInterior + ghost + 1; i++)
			{
				// Symmetry.
				u[IDX(i, k)] = (double)z_sym * u[IDX(i, 2 * ghost - 1 - k)];
			}
		}

		// Correct corner using diagonal symmetry.
		u[IDX(k, k)] = (double)(r_sym * z_sym) * u[IDX(2 * ghost - 1 - k, 2 * ghost - 1 - k)];
	}

	return;
}
