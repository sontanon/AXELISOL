// Allocate diff array.
void low_rank_allocate(const int ndiff);

// Deallocate diff array.
void low_rank_deallocate(void);

// Number of different elements in flat Laplacian matrices.
int ndiff_flat_laplacian(const int NrInterior, const int NzInterior);

// Number of different elements in general elliptic equation matrices.
int ndiff_general_elliptic(const int NrInterior, const int NzInterior, const int norder);

// Fill flat Laplacian diff array.
void low_rank_flat_laplacian(const int NrInterior, const int NzInterior);

// Fill general elliptic equation diff array.
void low_rank_general_elliptic(const int NrInterior, const int NzInterior, const int norder);
