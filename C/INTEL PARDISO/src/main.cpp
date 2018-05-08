// Global headers and variables.
#include "tools.h"

// PARDISO tools.
#include "pardiso_start.h"
#include "pardiso_stop.h"
#include "low_rank.h"

// Flat solver.
#include "flat_laplacian.h"

// General solver.
#include "general_elliptic.h"

// SOLVER RANGES.
#define NRINTERIOR_MIN 32
#define NRINTERIOR_MAX 2048
#define NZINTERIOR_MIN 32
#define NZINTERIOR_MAX 2048
#define DR_MAX 1.0
#define DR_MIN 0.000976562
#define DZ_MAX 1.0
#define DZ_MIN 0.000976562

int main(int argc, char *argv[])
{
	// PARAMETERS: Default values.
	int NrInterior = 256;
	int NzInterior = 64;
	double dr = 0.03125;
	double dz = 0.125;
	int norder = 2;
	char solver[256] = "flat";
	char dirname[256] = "output";
	int nrobin = 1;
	// Grid variables.
	int NrTotal = 0;
	int NzTotal = 0;
	int ghost = 0;
	int DIM = 0;
	// Various timers.
	clock_t start_time[10];
	clock_t end_time[10];
	double time[10];

	// User input character.
	char opt;

	// Get arguments from command line, first do a sanity check.
	if (argc < 8)
	{
		printf("ELLSOLVEC: WARNING! Usage is  $./ELLSOLVEC dirname solver norder NrInterior NzInterior dr dz\n");
		printf("           [solver] is the type of solver: flat or general.\n");
		printf("           [dirname] is a valid directory string name.\n");
		printf("           [norder] is an integer equal to 2 or 4 corresponding to the finite difference order.\n");
		printf("           [NrInterior] and [NzInterior] are integers equal to the number of interior points in r, z.\n");
		printf("           [dr] and [dz] are floating point doubles equal to the spatial step in r, z.\n");
		printf("           [nrobin] is an optional argument corresponding to Robin operator order: 1, 2, 3.\n");
		printf("Press (y/n) to procede with default arguments:\n");
		opt = getchar();
		getchar();
		if ((opt =='y') || (opt == 'Y'))
		{
			printf("ELLSOLVEC: User chose to proceed with default arguments.\n");
		}
		else
		{
			printf("ELLSOLVEC: User chose to abort.\n");
			exit(1);
		}
	}
	else
	{
		// Get arguments doing some sanity checks.
		// Directory name.
		memset(dirname, 0, 256);
		strcpy(dirname, argv[1]);

		// Solver type.
		memset(solver, 0, 256);
		strcpy(solver, argv[2]);
		if (!(strcmp(solver, "flat") == 0) && !(strcmp(solver, "general") == 0))
		{
			printf("ELLSOLVEC: ERROR! Unrecognized solver %s. Only \"flat\" or \"general\" is supported.\n", solver);
			exit(1);
		}

		// Finite difference order.
		norder = atoi(argv[3]);
		if ((norder != 2) && (norder != 4))
		{
			printf("ELLSOLVEC: ERROR! Finite difference %d is not supported, only 2 or 4.\n", norder);
			exit(1);
		}

		// Number of interior points.
		NrInterior = atoi(argv[4]);
		NzInterior = atoi(argv[5]);
		if ((NrInterior < NRINTERIOR_MIN) || (NzInterior > NRINTERIOR_MAX) || NrInterior <= 0)
		{
			printf("ELLSOLVEC: ERROR! NrInterior = %d is out of range [%d, %d].\n", NrInterior, NRINTERIOR_MIN, NRINTERIOR_MAX);
			exit(1);
		}
		if ((NzInterior < NZINTERIOR_MIN) || (NzInterior > NZINTERIOR_MAX) || NzInterior <= 0)
		{
			printf("ELLSOLVEC: ERROR! NzInterior = %d is out of range [%d, %d].\n", NzInterior, NZINTERIOR_MIN, NZINTERIOR_MAX);
			exit(1);
		}

		// Spatial steps.
		dr = atof(argv[6]);
		dz = atof(argv[7]);
		if ((dr < DR_MIN) || (dr > DR_MAX) || dr <= 0.0)
		{
			printf("ELLSOLVEC: ERROR! dr = %3.3E out of range [%3.3E, %3.3E].\n", dr, DR_MIN, DR_MAX);
			exit(1);
		}
		if ((dz < DZ_MIN) || (dz > DZ_MAX) || dz <= 0.0)
		{
			printf("ELLSOLVEC: ERROR! dz = %3.3E out of range [%3.3E, %3.3E].\n", dz, DZ_MIN, DZ_MAX);
			exit(1);
		}
	}
	// Get possible Robin operator.
	if (argc == 9)
	{
		nrobin = atoi(argv[8]);
		if (nrobin != 1 && nrobin != 2 && nrobin != 3)
		{
			printf("ELLSOLVEC: ERROR! nrobin = %d is not supported! Only 1, 2, 3.\n", nrobin);
			exit(1);
		}
	}

	// Do I/O on output directory.
	// First create directory.
	struct stat st = { 0 };
	if (stat(dirname, &st) == -1)
	{
#ifdef WIN
		_mkdir(dirname);
#else
		mkdir(dirname, 0755);
#endif
	}
	else 
	{
		printf("ELLSOLVEC: WARNING! Directory %s already exists.\n", dirname);
		printf("Press (y/n) to procede and possibly overwrite files:\n");
		opt = getchar();
		getchar();
		if ((opt =='y') || (opt == 'Y'))
		{
			printf("ELLSOLVEC: User chose to proceed.\n");
		}
		else
		{
			printf("ELLSOLVEC: User chose to abort.\n");
			exit(1);
		}
	}
	// Now CD to output directory.
	if (chdir(dirname) == -1)
	{
		printf("ELLSOLVEC: WARNING! Could not CD to %s directory.\n", dirname);
		printf("Press (y/n) to procede and write in current directory:\n");
		opt = getchar();
		getchar();
		if ((opt =='y') || (opt == 'Y'))
		{
			printf("ELLSOLVEC: User chose to proceed.\n");
		}
		else
		{
			printf("ELLSOLVEC: User chose to abort.\n");
			exit(1);
		}
	}

	// First get finite difference order.
	if (norder == 2)
	{
		ghost = 2;
	}
	else if (norder == 4)
	{
		ghost = 3;
	}
	// Calculate longitudinal dimensions.
	NrTotal = ghost + NrInterior + 1;
	NzTotal = ghost + NzInterior + 1;
	DIM = NrTotal * NzTotal;

	// Print info to screen.
	printf("ELLSOLVEC: System parameters are:\n");
	printf("\tdirname\t= %s\n", dirname);
	printf("\tsolver\t= %s\n", solver);
	printf("\torder\t= %d\n", norder);
	printf("\tNrInterior\t= %d\n", NrInterior);
	printf("\tNzInterior\t= %d\n", NzInterior);
	printf("\tdr\t= %4.8E\n", dr);
	printf("\tdz\t= %4.8E\n", dz);
	printf("\tnrobin\t= %d\n", nrobin);

	// Grid functions.
	double *r, *z, *u, *f, *s, *res;
	double *a, *b, *c, *d, *e;
	// Size of double memory.
	size_t DIM_size = DIM * sizeof(double);

	// Allocate variables.
	printf("ELLSOLVEC: Allocating memory...\n");
	r = (double *)malloc(DIM_size);
	z = (double *)malloc(DIM_size);
	u = (double *)malloc(DIM_size);
	f = (double *)malloc(DIM_size);
	s = (double *)malloc(DIM_size);
	res = (double *)malloc(DIM_size);
	a = (double *)malloc(DIM_size);
	b = (double *)malloc(DIM_size);
	c = (double *)malloc(DIM_size);
	d = (double *)malloc(DIM_size);
	e = (double *)malloc(DIM_size);
	printf("ELLSOLVEC: Allocated memory.\n");
	
	// Auxiliary variables.
	double aux_r, aux_z, aux1, aux2;
	int i, j, k;

	// Fill grids.
	#pragma omp parallel shared(r, z, u, f, s, res, a, b, c, d, e) private(aux_r, aux_z, j, k)
	{
		#pragma omp for schedule(guided)
		for (i = 0; i < NrTotal; i++)
		{
			aux_r = ((double)(i - ghost) + 0.5) * dr;
			for (j = 0; j < NzTotal; j++)
			{
				aux_z = ((double)(j - ghost) + 0.5) * dz;
				k = IDX(i, j);
				r[k] = aux_r;
				z[k] = aux_z;
				// Set everything else to zero.
				u[k] = 0.0;
				f[k] = 0.0;
				s[k] = 0.0;
				res[k] = 0.0;
				a[k] = 0.0;
				b[k] = 0.0;
				c[k] = 0.0;
				d[k] = 0.0;
				e[k] = 0.0;
			}
		}
	}
	// Write coordinate grids.
	write_single_file(r, "r.asc", NrTotal, NzTotal);
	write_single_file(z, "z.asc", NrTotal, NzTotal);
	printf("ELLSOLVEC: Filled and printed r, z grids.\n");

	// Intialize memory and parameters.
	pardiso_start(NrInterior, NzInterior);

	// Choose between flat and general solver.
    	// Flat solver.
	if (strcmp(solver, "flat") == 0)
	{
		// Fill linear source and RHS.
		#pragma omp parallel shared(s, f) private(aux_r, aux_z)
		{
			#pragma omp for schedule(guided)
			for (k = 0; k < DIM; k++)
			{
				// Coordinates.
				aux_r = r[k];
				aux_z = z[k];
				// Linear source.
				s[k] = exp(-aux_r * aux_r - aux_z * aux_z) * (0.5 + aux_r * aux_r * (-3.0 + aux_r * aux_r + aux_z * aux_z));
				// RHS.
				f[k] = 0.0;
			}
		}

		// Write linear source and RHS.
		write_single_file(s, "s.asc", NrTotal, NzTotal);
		write_single_file(f, "f.asc", NrTotal, NzTotal);

		// Call solver.
		printf("ELLSOLVEC: Calling normal solver.\n");
		start_time[0] = clock();
		flat_laplacian(u, res, s, f, 1.0, nrobin, 1, 1, 
			NrInterior, NzInterior, ghost, dr, dz, norder,
			0, 0);
		end_time[0] = clock();
		time[0] = (double)(end_time[0] - start_time[0])/CLOCKS_PER_SEC;

		// Precondition with CGS.
		printf("ELLSOLVEC: Solving with CGS.\n");
		start_time[3] = clock();
		flat_laplacian(u, res, s, f, 1.0, nrobin, 1, 1, 
			NrInterior, NzInterior, ghost, dr, dz, norder,
			0, 6);
		end_time[3] = clock();
		time[3] = (double)(end_time[3] - start_time[3])/CLOCKS_PER_SEC;

		// Low rank update solve.
		// Get number of differing elements.
		int ndiff = ndiff_flat_laplacian(NrInterior, NzInterior);
		printf("ELLSOLVEC: Number of differing elements in low rank update are %d.\n", ndiff);
		// Allocate diff array.
		low_rank_allocate(ndiff);
		// Fill diff array for flat Laplacian.
		low_rank_flat_laplacian(NrInterior, NzInterior);

		// Call solver with low rank update.
		printf("ELLSOLVEC: Solving with low rank update.\n");
		start_time[5] = clock();
		flat_laplacian(u, res, s, f, 1.0, nrobin, 1, 1, 
			NrInterior, NzInterior, ghost, dr, dz, norder,
			1, 0);
		end_time[5] = clock();
		time[5] = (double)(end_time[5] - start_time[5])/CLOCKS_PER_SEC;
	}
    	// General solver.
	else if (strcmp(solver, "general") == 0)
	{
		// Fill coefficients, linear source and RHS.
		#pragma omp parallel shared(a, b, c, d, e, s, f) private(aux_r, aux_z)
		{
			#pragma omp for schedule(guided)
			for (k = 0; k < DIM; k++)
			{
				// Coordinates.
				aux_r = r[k];
				aux_z = z[k];
				// Coefficients.
				a[k] = aux_r;
				b[k] = 0.0;
				c[k] = aux_r;
				d[k] = 1.0;
				e[k] = 0.0;
				// Linear source.
				s[k] = aux_r * exp(-aux_r * aux_r - aux_z * aux_z) * (0.5 + aux_r * aux_r * (-3.0 + aux_r * aux_r + aux_z * aux_z));
				// RHS.
				f[k] = 0.0;
			}
		}
		// Write coefficients.
		write_single_file(a, "a.asc", NrTotal, NzTotal);
		write_single_file(b, "b.asc", NrTotal, NzTotal);
		write_single_file(c, "c.asc", NrTotal, NzTotal);
		write_single_file(d, "d.asc", NrTotal, NzTotal);
		write_single_file(e, "e.asc", NrTotal, NzTotal);
		write_single_file(s, "s.asc", NrTotal, NzTotal);
		write_single_file(f, "f.asc", NrTotal, NzTotal);

		// Call solver.
		printf("ELLSOLVEC: Calling normal solver.\n");
		start_time[0] = clock();
		general_elliptic(u, res, a, b, c, d, e, s, f, 1.0, nrobin, 1, 1, 
			NrInterior, NzInterior, ghost, dr, dz, norder,
			0, 0);
		end_time[0] = clock();
		time[0] = (double)(end_time[0] - start_time[0])/CLOCKS_PER_SEC;

		// Precondition with CGS.
		printf("ELLSOLVEC: Solving with CGS.\n");
		start_time[3] = clock();
		general_elliptic(u, res, a, b, c, d, e, s, f, 1.0, nrobin, 1, 1, 
			NrInterior, NzInterior, ghost, dr, dz, norder,
			0, 6);
		end_time[3] = clock();
		time[3] = (double)(end_time[3] - start_time[3])/CLOCKS_PER_SEC;

		// Low rank update solve.
		// Get number of differing elements.
		int ndiff = ndiff_general_elliptic(NrInterior, NzInterior, norder);
		printf("ELLSOLVEC: Number of differing elements in low rank update are %d.\n", ndiff);
		// Allocate diff array.
		low_rank_allocate(ndiff);
		// Fill diff array for general elliptic.
		low_rank_general_elliptic(NrInterior, NzInterior, norder);

		// Call solver with low rank update.
		printf("ELLSOLVEC: Solving whith low rank update.\n"); 
		start_time[5] = clock(); 
		general_elliptic(u, res, a, b, c, d, e, s, f, 1.0, nrobin, 1, 1, 
		NrInterior, NzInterior, ghost, dr, dz, norder,
			1, 0);
		end_time[5] = clock();
		time[5] = (double)(end_time[5] - start_time[5])/CLOCKS_PER_SEC;
	}

	// Print execution times.
	printf("ELLSOLVEC: Normal solver took %3.3E seconds.\n", time[0]);
	printf("ELLSOLVEC: Solver with CGS took %3.3E seconds.\n", time[3]);
	printf("ELLSOLVEC: Solver with low rank update took %3.3E seconds.\n", time[5]);

	// Write solution and residual.
	write_single_file(u, "u.asc", NrTotal, NzTotal);
	write_single_file(res, "res.asc", NrTotal, NzTotal);

	// Deallocate low rank array: same for flat or general solver.
	low_rank_deallocate();

	// Clear memory and parameters.
	pardiso_stop();

	// Deallocate memory.
	printf("ELLSOLVEC: Cleaning up...\n");
	free(r);
	free(z);
	free(u);
	free(f);
	free(s);
	free(res);
	free(a);
	free(b);
	free(c);
	free(d);
	free(e);
	printf("ELLSOLVEC: Cleared all memory.\n");

	// All done.
	return 0;
}
