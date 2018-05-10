// Allocate diff array.
void low_rank_allocate(const int ndiff);

// Deallocate diff array.
void low_rank_deallocate(void);

// Number of different elements between matrices.
// Flat Laplacian.
int ndiff_flat_laplacian(const int NrInterior, const int NzInterior);

// General elliptic equation.
int ndiff_general_elliptic(const int NrInterior, const int NzInterior, const int order);

// Fill diff array.
// Flat Laplacian.
void low_rank_flat_laplacian(const int NrInterior, const int NzInterior);

// General elliptic equation.
void low_rank_general_elliptic(const int NrInterior, const int NzInterior, const int order);
