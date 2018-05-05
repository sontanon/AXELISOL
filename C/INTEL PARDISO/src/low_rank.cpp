// One-based indexing is in tools header.
#include "tools.h"
#include "pardiso_param.h"


// Number of different elements between matrices.
// Flat Laplacian.
int ndiff_flat_laplacian(const int NrInterior, const int NzInterior)
{
	// Number of elements is independent of order.
	int ndiff = NrInterior * NzInterior;
	
	return ndiff;
}

// Allocate diff array.
void low_rank_allocate(const int ndiff)
{
	// Sanity check that diff array has not already been allocated.
	if (diff != NULL)
	{
		printf("ERROR: Low Rank diff array was already allocated!\n");
		exit(1);
	}

	// Allocate memory for diff array.
	diff = (int *)malloc((2 * ndiff + 1) * sizeof(int));

#ifdef VERBOSE
	printf("PARDISO LOW RANK UPDATE: Setup diff and ndiff parameters.\n");
#endif
}

// Deallocate diff array.
void low_rank_deallocate(void)
{
	// Deallocate.
	free(diff);

	// Point towards NULL.
	diff = NULL;
}

// Fill diff array.
// Flat Laplacian.
void low_rank_flat_laplacian(const int NrInterior, const int NzInterior)
{
	// Auxiliary integers.
	int i, j;

	// Number of elements we have filled in.
	int offset = 0;

	// Temporary offset.
	int t_offset = 0;

	// First element of diff array is ndiff itself.
	diff[offset] = ndiff_flat_laplacian(NrInterior, NzInterior);

	// Increase offset.
	offset += 1;

	// Set temporary offset.
	t_offset = offset;

	// Different elements are the same independent of order.
	// They are always the interior diagonal elements.
	#pragma omp parallel shared(diff) private(offset, j)
	{
		#pragma omp for schedule(guided)
		for (i = 1; i < NrInterior + 1; i++)
		{
			// Each iteration of i loop will fill 2 * NzInterior elements in diff array. 
			offset = t_offset + (i - 1) * 2 * NzInterior;

			// Loop over interior points.
			for (j = 1; j < NzInterior + 1; j++)
			{
				// Row indices.
				diff[offset] = BASE + IDX(i, j);
				// Column indices.
				diff[offset + 1] = BASE + IDX(i, j);
				// Increase offset.
				offset += 2;
			}
		}
	}

	// All done.
	return;
}
