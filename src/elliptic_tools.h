// Reduce array u to the elliptic solver-sized g_u.
void ghost_reduce(const double *u, double *g_u, const int NrInterior, const int NzInterior, const int ghost);

// Fill array u from elliptic solver sized-array g_u using symmetry conditions.
void ghost_fill(const double *g_u, double *u, const int r_sym, const int z_sym, const int NrInterior, const int NzInterior, const int ghost);
