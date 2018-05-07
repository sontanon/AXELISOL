// Allocate diff array.
extern "C" void low_rank_allocate(const int ndiff);

// Deallocate diff array.
extern "C" void low_rank_deallocate(void);

// Number of different elements in flat Laplacian matrices.
extern "C" int ndiff_flat_laplacian(const int NrInterior, const int NzInterior);

// Number of different elements in general elliptic equation matrices.
extern "C" int ndiff_general_elliptic(const int NrInterior, const int NzInterior, const int norder);

// Fill flat Laplacian diff array.
extern "C" void low_rank_flat_laplacian(const int NrInterior, const int NzInterior);

// Fill general elliptic equation diff array.
extern "C" void low_rank_general_elliptic(const int NrInterior, const int NzInterior, const int norder);
