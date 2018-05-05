// Number of different elements in flat Laplacian matrices.
int ndiff_flat_laplacian(const int NrInterior, const int NzInterior);

// Allocate diff array.
void low_rank_allocate(const int ndiff);

// Deallocate diff array.
void low_rank_deallocate(void);

// Fill flat Laplacian diff array.
void low_rank_flat_laplacian(const int NrInterior, const int NzInterior);
