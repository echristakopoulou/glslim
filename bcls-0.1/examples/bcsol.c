// =====================================================================
// $Revision: 290 $ $Date: 2007-03-04 22:01:54 -0800 (Sun, 04 Mar 2007) $
//
// Driver for BCLS: a bound-constrained least-squares solver.
//
// 25 Aug 05: Original version.
//            Michael P. Friedlander
//            mpf@cs.ubc.ca
//            University of British Columbia
// =====================================================================

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <getopt.h>
#include "cs.h"    // The CSparse matrix library.
#include "iohb.h"  // The Harwell-Boing file reader library.
#include "bcls.h"  // The BCLS solver library.

// ---------------------------------------------------------------------
// Global variables.
// ---------------------------------------------------------------------

// Name of the (required) input file for A and b.
static char *in_file_Ab_name = NULL;

// Name of the (optional) input file for bl, bu, c, x.
static char *in_file_opt_name = NULL;

// Name of the (optional) minor-output file.
static char *out_file_minor = NULL;

// Optional storage description.  At least RHS must be stored.
static char opt_storage[5] = "lu";

// Minor iteration limit.
static int minor_its = -1;

// Major iteration limit.
static int major_its = -1;

// Optimality tolerance.
static double opt_tol = -1.0;

// Preconditioning.  None by default.
static int preconditioning = 0;

// Print level.
static int print_level = -1;

// Type of projected search.
static int proj_search = BCLS_PROJ_SEARCH_EXACT;

// Method for Newton step computation.
static int newton_step = BCLS_NEWTON_STEP_LSQR;

// Method for Newton step computation.
static int scaled_steepest = 0;

// Regularization parameter.
static double damp = 0.0;

// If this flag is set, only check input, don't solve the problem.
static int check = 0;

// ---------------------------------------------------------------------
// Workspace structure passed to Aprod and Usolve routines.
// ---------------------------------------------------------------------
typedef struct {
    cs *A;
} worksp;

// ---------------------------------------------------------------------
// xmalloc
// Malloc wrapper.
//
// Returns a pointer to the newly allocated memory, or NULL if malloc
// returned an eror.
// ---------------------------------------------------------------------
static void *
xmalloc(int n, size_t size, char * who) {

    register void *value = malloc(n * size);
    if (value == NULL) {
	fprintf( stderr, "Not enough memory to allocate %s.\n", who );
	exit(EXIT_FAILURE);
    }
    return value;
}

// ---------------------------------------------------------------------
// dload
// Load a constant into a vector.
// ---------------------------------------------------------------------
static void
dload( const int n, const double alpha, double x[] ) {

    int j;
    for (j = 0; j < n; j++) x[j] = alpha;
    return;
}

// ---------------------------------------------------------------------
// display_help
// Display help at the commandline.
// ---------------------------------------------------------------------
static void
display_help(char *my_name) {

    printf("Usage: %s [options...] filename\n", my_name);
    printf("The last argument is the filename of a Harwell-Boeing file\n"
	   "that contains A and b (1-based indexing).\n");
    printf("Options:\n");
    printf("   -c, --check        do not solve problem, check input"
	   " data only\n");
    printf("       --cgls         use CGLS for Newton step computation\n");
    printf("   -d, --damp         damping (regularization) parameter\n");
    printf("   -a, --approx       approximate projected linesearch\n");
    printf("   -h, --help         print this message\n");
    printf("   -m, --minor        minor (CGLS/LSQR) iteration limit\n");
    printf("   -M, --major        major iteration limit\n");
    printf("       --minor-out    minor (CGLS/LSQR) output file."
           "Use - for stdout.\n");
    printf("   -O  --optional     optional data file for l, u, c, x\n");
    printf("   -o  --opttol       optimality tolerance\n");
    printf("       --precon       use LU preconditioner\n");
    printf("   -p, --print        print level 0-6 (default is 1)\n");
    printf("       --scaled       scaled steepest descent"
           " (default is unscaled)\n");
    printf("   -s, --storage      storage description in optional data file"
           " (default is 'lu')\n");
    return;
}

// ---------------------------------------------------------------------
// parse_cmdline
// Parse the commandline options.
// ---------------------------------------------------------------------
static void
parse_cmdline(int argc, char *argv[]) {

    int c, j;
    int option_index;
    static struct option long_options[] = {
        // Long options w/o short equivalents.  Only set a flag.
        {"cgls",     no_argument, &newton_step,   BCLS_NEWTON_STEP_CGLS},
        {"scaled",   no_argument, &scaled_steepest, 1},
        {"precon",   no_argument, &preconditioning, 1},
        // Long and short options.
        {"check",    no_argument,       0, 'c'},
        {"damp",     required_argument, 0, 'd'},
        {"approx",   no_argument,       0, 'a'},
        {"help",     no_argument,       0, 'h'},
        {"minor",    required_argument, 0, 'm'},
        {"major",    required_argument, 0, 'M'},
        {"minor-out",required_argument, 0, 'z'},
        {"optional", required_argument, 0, 'O'},
        {"opttol",   required_argument, 0, 'o'},
        {"print",    required_argument, 0, 'p'},
        {"storage",  required_argument, 0, 's'},
        // End of option flags.
        {0,0,0,0}
    };
    
    while (1) {
        
        c = getopt_long(argc, argv, "cd:ahm:M:O:o:p:s:",
                        long_options, &option_index );
        
        if ( c == -1 ) // No more options.
            break;

        switch (c) {

        case 0:   // Long option with no argument.
            break;
            
        case 'c': // -c  --check
            check = 1;
            break;

        case 'd': // -d  --damp
            damp = atof(optarg);
            break;
            
        case 'a': // -a --approx
            proj_search = BCLS_PROJ_SEARCH_APPROX;
            break;

        case 'h': // -h  --help
            display_help(argv[0]);
            exit(EXIT_SUCCESS);
            break;

        case 'm': // -m  --minor
            minor_its = atoi(optarg);
            break;

        case 'M': // -M  --major
            major_its = atoi(optarg);
            break;

        case 'O': // -O  --optional
            if (in_file_opt_name != NULL) {
                printf("Only one optional data file allowed\n");
                exit(EXIT_FAILURE);
            }
            in_file_opt_name = optarg;
            break;

        case 'z': // --minor-out
            out_file_minor = optarg;
            break;

        case 'o': //     --opttol
            opt_tol = atof(optarg);
            break;

        case 'p': // -p  --print
            print_level = atoi(optarg);
            if (print_level < 0) {
                printf("print level must be 0 or positive");
                exit(EXIT_FAILURE);
            }
            break;

        case 's': // -s  --storage
	    strncpy( opt_storage, optarg, sizeof(opt_storage) );
	    for (j = 0; j < strlen(opt_storage); j++) {
		if (opt_storage[j] != 'l' &&
		    opt_storage[j] != 'u' &&
		    opt_storage[j] != 'c' &&
		    opt_storage[j] != 'x' ) {
		    printf("Storage descriptor '%s' not recognized.\n",
			   &(opt_storage[j]));
		    exit(EXIT_FAILURE);
		}
	    }
            break;

        case '?':
            // getopt_long already printed an error message.
            printf("Try %s --help\n", argv[0]);
            exit(EXIT_FAILURE);
            break;

        default:
            display_help(argv[0]);
            exit(EXIT_SUCCESS);
            
        }
    }
    
    // Gather the remaining non-option arguments.
    if (optind < argc) {
        while (optind < argc) {
            if (in_file_Ab_name != NULL) {
                printf("Only one data input file allowed\n");
                printf("Try %s --help\n", argv[0]);
                exit(EXIT_FAILURE);
            }
            in_file_Ab_name = argv[optind];
            optind++;
        }
    }
    else {
        printf("No input file for A, b specified; try %s --help\n", argv[0]);
        exit(EXIT_FAILURE);
    }
}

// ---------------------------------------------------------------------
// Aprod
// Matrix-vector products.
//
// If     mode == BCLS_PROD_A  (0), compute y <- A *x, with x untouched;
// and if mode == BCLS_PROD_At (1), compute x <- A'*y, with y untouched.
// ---------------------------------------------------------------------
static int
Aprod( int mode, int m, int n, int nix, int ix[],
       double x[], double y[], void *UsrWrk ) {

    int     i, j, k, l;
    double  aij, xj, sum;
    worksp * Wrk = (worksp *)UsrWrk;
    cs *A = (cs *)Wrk->A;

    assert( mode == BCLS_PROD_A     ||
            mode == BCLS_PROD_At    ||
            mode == BCLS_PROD_INIT  ||
            mode == BCLS_PROD_TERM  );
    assert( nix  <= n );
    assert( A->m == m );
    assert( A->n == n );

    if (mode == BCLS_PROD_A) {

        dload( m, 0.0, y );

	for (l = 0; l < nix; l++) {
	    j = ix[l];
	    xj = x[j];
	    if (xj == 0.0)
		; // Relax.
	    else
		for (k = A->p[j]; k < A->p[j+1]; k++) {
		    aij   = A->x[k];
		    i     = A->i[k];
		    y[i] += aij * xj;
		}
	}
    }

    else if (mode == BCLS_PROD_At) {
	for (l = 0; l < nix; l++) {
	    j = ix[l];
	    sum = 0;
	    for (k = A->p[j]; k < A->p[j+1]; k++) {
		aij  = A->x[k];
		i    = A->i[k];
		sum += aij * y[i];
	    }
	    x[j] = sum;
	}
    }

    // Exit.
    return 0;
}

// ---------------------------------------------------------------------
// Usolve
// Preconditioner.
//
// If     mode = BCLS_PRECON_U,  solve  U v = w,
// and if mode = BCLS_PRECON_Ut, solve  U'w = v.
// ---------------------------------------------------------------------
static int
Usolve( int mode, int m, int n, int nix, int ix[],
        double v[], double w[], void *UsrWrk ){
    
    assert( nix  <= n );
    assert( mode == BCLS_PRECON_U    ||
            mode == BCLS_PRECON_Ut   ||
            mode == BCLS_PRECON_INIT ||
            mode == BCLS_PRECON_TERM );
    
    if ( mode == BCLS_PRECON_INIT ) {
        // -------------------------------------------------------------
        // Initialize the preconditioner.
        // -------------------------------------------------------------

        // -------------------------------------------------------------
        // Factorize A(:,ix).
        // -------------------------------------------------------------


    }
    else if ( mode == BCLS_PRECON_U ) {
        
        // Solve  U v = w.

    }
    else if ( mode == BCLS_PRECON_Ut ) {

        // Solve  U'w = v.

    }
    else if ( mode == BCLS_PRECON_TERM )
        ; // Relax.
    
    // Exit.
    return 0;
}

// ---------------------------------------------------------------------
// CallBack.
// Periodically called by BCLS to test if the user wants to exit.
// ---------------------------------------------------------------------
static int
CallBack( BCLS *ls, void *UsrWrk )
{
    int err;
    err = 0;  // No error.
    return err;
}

// ---------------------------------------------------------------------
// pretty_printer
// This is the print-routine that will be used by BCLS for its output.
// ---------------------------------------------------------------------
static int
pretty_printer( void *io_file, char *msg ) {
    fprintf( io_file, msg );
    return 0;
}

// ---------------------------------------------------------------------
// main
// ---------------------------------------------------------------------
int
main(int argc, char *argv[]) {

    // Miscellaneous variables.
    int i, j, err, nRhs;
    double *dPtr;
    char line[1025], *blActiv, *buActiv, *MatType;
    FILE *in_file;

    // These variables help to define the BCLS problem
    BCLS   *ls;              // A BCLS problem.
    worksp  Wrk;             // Workspace.
    int     m, n, nnz;       // Problem dimensions.
    double *c     = NULL;    // Linear vector.
    double *x     = NULL;    // Solution vector.
    double *b     = NULL;    // RHS vector.
    double *bl    = NULL;    // Lower bounds.
    double *bu    = NULL;    // Upper bounds.
    double *anorm = NULL;    // Column norms of A.

    // -----------------------------------------------------------------
    // Basic initialization.
    // -----------------------------------------------------------------
    parse_cmdline( argc, argv );

    // -----------------------------------------------------------------
    // Read in the header for A and b.
    // -----------------------------------------------------------------
    if (!readHB_info( in_file_Ab_name, &m, &n, &nnz, &MatType, &nRhs ))
        exit(EXIT_FAILURE);
    free(MatType);
    printf( " ----------------------------------------------------\n");
    printf( " Reading from data file %s:\n", in_file_Ab_name );
    printf( " Matrix size:  m = %i  n = %i  nnz = %i  nRhs = %i\n",
            m, n, nnz, nRhs );
    printf( " ----------------------------------------------------\n");

    // -----------------------------------------------------------------
    // Allocate storage for A, b, bl, bu, x, c.
    // -----------------------------------------------------------------

    // Allocate storage for A.
    Wrk.A = cs_spalloc( m, n, nnz, 1, 0 );

    // Storage for these vectors are always needed.
    b  = (double *)xmalloc( m, sizeof(double), "b"  );
    bl = (double *)xmalloc( n, sizeof(double), "bl" );
    bu = (double *)xmalloc( n, sizeof(double), "bu" );
    x  = (double *)xmalloc( n, sizeof(double), "x"  );

    // Storage for c is optional.
    if (strchr(opt_storage, 'c'))
    c  = (double *)xmalloc( n, sizeof(double), "c"  );

    // -----------------------------------------------------------------
    // Load A, b.
    // -----------------------------------------------------------------
    readHB_mat_double( in_file_Ab_name, Wrk.A->p, Wrk.A->i, Wrk.A->x );
    readHB_aux_double( in_file_Ab_name, 'F', b );

    // -----------------------------------------------------------------
    // Load bl, bu, c, x.
    // -----------------------------------------------------------------

    // Load bl, bu, x with defaults if no values are specified.
    if (!strchr(opt_storage,'l') || !in_file_opt_name) dload(n,-1e20,bl);
    if (!strchr(opt_storage,'u') || !in_file_opt_name) dload(n, 1e20,bu);
    if (!strchr(opt_storage,'x') || !in_file_opt_name) dload(n,  0.0, x);

    // Open optional file for reading.
    if (in_file_opt_name) {
        if ( (in_file = fopen( in_file_opt_name, "r")) == NULL) {
            fprintf(stderr,"Error: Cannot open file: %s\n",
                    in_file_opt_name);
            exit(EXIT_FAILURE);
        }
        
        for (j = 0; j < strlen(opt_storage); j++) {
            
            if      (opt_storage[j] == 'l')  dPtr = bl;
            else if (opt_storage[j] == 'u')  dPtr = bu;
            else if (opt_storage[j] == 'x')  dPtr = x;  
            else                             dPtr = c;
            
            for (i = 0; i < n; i++ ) {
                if ( fgets( line, sizeof(line), in_file ) == NULL ) {
                    fprintf( stderr, "Reached end of file %s before "
                             "reading %d entries in mode %s.\n",
                             in_file_opt_name, i, &(opt_storage[j]) );
                    exit(EXIT_FAILURE);
                }
                sscanf( line, "%lg", &(dPtr[i]) );
            }
        }
        
        // Done with input file.
        fclose( in_file );
    }
    // -----------------------------------------------------------------
    // Initialize a BCLS problem.  This routine MUST be called before
    // any other BCLS routine.
    // -----------------------------------------------------------------
    ls = bcls_create_prob( m, n );

    // Compute column norms.
    if ( scaled_steepest ) {

        // Allocate memory to hold the scales.
        anorm = xmalloc( n, sizeof(double), "anorm" );

        // Compute the scales and let BCLS know we have them.
        bcls_compute_anorm( ls, n, m, Aprod, &Wrk, anorm );
        bcls_set_anorm( ls, anorm );
    }

    // Instatiate a particular BCLS problem.
    bcls_set_problem_data( ls,    // The BCLS problem
			   m,     // Number of problem rows
			   n,     // Number of problem columns
			   Aprod, // The Mat-vec routine
			   &Wrk,  // Arbitrary data for the Mat-vec routine
			   damp,  // Damping parameter
			   x,     // Solution vector
			   b,     // RHS vector
			   c,     // Linear term (may be NULL)
			   bl,    // Lower-bounds vector
			   bu );  // Upper-bounds vector

    // Set the user options.
    bcls_set_print_hook( ls, stdout, pretty_printer );
    if (major_its   > 0) ls->itnMajLim     = major_its;
    if (minor_its   > 0) ls->itnMinLim     = minor_its;
    if (opt_tol     > 0) ls->optTol        = opt_tol;
    if (print_level > 0) ls->print_level   = print_level;
    if (preconditioning) bcls_set_usolve( ls, Usolve );
    if (out_file_minor) {
        if (strcmp(out_file_minor, "-") == 0)
            ls->minor_file = stdout;
        else
            ls->minor_file = fopen(out_file_minor, "w+");
    }
    ls->proj_search    = proj_search;
    ls->newton_step    = newton_step;
    ls->CallBack       = CallBack;
                         
    // Call the main BCLS routine.
    err = bcls_solve_prob( ls );

    // Print the solution if print level >= 4.
    if (print_level >= 4) {
	printf("\n Solution\n --------\n");
	printf("%4s  %18s %1s %18s %1s %18s  %18s\n",
	       "Var","Lower","","Value","","Upper","Gradient");
	for (j = 0; j < n; j++) {
	    blActiv = "";
	    buActiv = "";
	    if (x[j] - bl[j] < ls->epsx) blActiv = "=";
	    if (bu[j] - x[j] < ls->epsx) buActiv = "=";
 	    printf("%4d  %18.11e %1s %18.11e %1s %18.11e  %18.11e\n",
		   j+1, bl[j], blActiv, x[j], buActiv, bu[j], (ls->g)[j]);
	}
    }
    // Deallocate the BCLS problem.
    err = bcls_free_prob( ls );

    // Deallocate memory.
    cs_spfree( Wrk.A );
    if ( b   ) free( b     );
    if ( bl  ) free( bl    );
    if ( bu  ) free( bu    );
    if ( x   ) free( x     );
    if ( c   ) free( c     );
    if (anorm) free( anorm );

    // Close files.
    if (out_file_minor)
        fclose( ls->minor_file );

    // Exit with no errors.  Whew!
    return (EXIT_SUCCESS);
}
