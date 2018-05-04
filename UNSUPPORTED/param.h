#ifdef MAIN_FILE
/* GRID */
double dr = 0.03125;
double dz = 0.125;
int NrInterior = 256;
int NzInterior = 64;
int NrTotal = 0;
int NzTotal = 0;
int ghost = 0;
int DIM = 0;
const char *order = "four";

/*
// OUTPUT //
const char *dirname = "Test";

// ELLIPTIC SOLVER //
int use_permutation = 0;
int use_preconditioner = 0;
int perm_flag = 0;
int perm_count = 0;
int perm_count_max = 0;
int precond_flag = 0;
int precond_L = 0;
int nrobin = 1;
*/
#else
/* GRID */
extern double dr;
extern double dz;
extern int NrInterior;
extern int NzInterior;
extern int NrTotal;
extern int NzTotal;
extern int ghost;
extern int DIM;
extern const char *order;

/*
// OUTPUT //
extern const char *dirname;

// ELLIPTIC SOLVER //
extern int use_permutation;
extern int use_preconditioner;
extern int perm_flag;
extern int perm_count;
extern int perm_count_max;
extern int precond_flag;
extern int precond_L;
extern int nrobin;
*/
#endif
