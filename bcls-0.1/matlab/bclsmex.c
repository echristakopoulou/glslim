/*======================================================================
   $Revision: 284 $ $Date: 2006-12-17 20:55:25 -0800 (Sun, 17 Dec 2006) $
   Matlab MEX interface to BCLS solver.
 
   13 Oct 05: Original version.
              Michael P. Friedlande
              mpf@cs.ubc.ca
              University of British Columbia
=======================================================================*/

#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mex.h"
#include "bcls.h"
#include "cblas.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

/* Flag if control-c was pressed. */
volatile sig_atomic_t ctrl_c_pressed = 0;
void (*MATLAB_handler)(int);
void (*BCLS_handler)(int);

/* =====================================================================
   Insert revision info into the object file.
   =====================================================================*/
static char *revisionid =
    "$Revision: 284 $"
    "$Date: 2006-12-17 20:55:25 -0800 (Sun, 17 Dec 2006) $";

/* =====================================================================
   Matlab helper routines.
   =====================================================================*/
typedef struct {  /* Matlab full matrix. */
  int m, n;
  double *val;
} matlab_Mat;

typedef struct {  /* Compressed column sparse format. */
  int    n, m, nnz;
  int    *ind, *col;
  double *val;
} matlab_spMat;

static double
assertScalar( const mxArray *mex, char who[8] )
{
   char msgbuf[60];
   int  m = mxGetM( mex );
   int  n = mxGetN( mex );
   if ( m != 1 || n != 1 || mxIsSparse( mex ) ) {
      sprintf( msgbuf,
	       "Expected '%.8s' to be a scalar.  Instead it's %dx%d.\n",
	       who,m,n);
      mexErrMsgTxt( msgbuf );
   }
   else {
      return mxGetScalar( mex );
   }
}

static char *
assertString( const mxArray *mex, char who[8], int add )
{
  char msgbuf[60];
  int  m = mxGetM( mex );
  int  n = mxGetN( mex );
  if ( m != 1 || mxIsSparse( mex ) ) {
    sprintf( msgbuf,
	     "Expected '%.8s' to be a single string.  Instead it's %dx%d.\n",
	     who,m,n);
    mexErrMsgTxt( msgbuf );
  }
  else {
    int  buflen  = mxGetN( mex ) + 1;
    char *string = mxCalloc( buflen + add, sizeof(char) );
    int err      = mxGetString( mex, string, buflen );
    if (err) {
      sprintf( msgbuf, "Not enough space.  '%.8s' was trucated.\n", who );
      mexWarnMsgTxt( msgbuf );
    }
    else {
      return string;
    }
  }
}

static double *
assertSpMatrix( const mxArray *mex,
		int m, int n, char who[8] )
{
   char msgbuf[60];
   int  m_actual = mxGetM( mex );
   int  n_actual = mxGetN( mex );
   if ( m != m_actual || n != n_actual || !mxIsSparse( mex ) ) {
      sprintf( msgbuf,
       "Expected '%.8s' to be a sparse %dx%d matrix.  Instead it's %dx%d.\n",
       who, m,n, m_actual,n_actual);
      mexErrMsgTxt( msgbuf );
   }
   else {
      return mxGetPr( mex );
   }
}

static double *
assertMatrix( const mxArray *mex,
	      int m, int n, char who[8] )
{
   char msgbuf[60];
   int  m_actual = mxGetM( mex );
   int  n_actual = mxGetN( mex );
   if ( m != m_actual || n != n_actual || mxIsSparse( mex ) ) {
      sprintf( msgbuf,
       "Expected '%.8s' to be a dense %dx%d matrix.  Instead it's %dx%d.\n",
       who, m,n, m_actual,n_actual);
      mexErrMsgTxt( msgbuf );
   }
   else {
      return mxGetPr( mex );
   }
}

static matlab_spMat
create_matlab_spmat( mxArray **plhs, int m, int n, int nnz )
{
  matlab_spMat A;
  A.n   = n;
  A.m   = m;
  A.nnz = nnz;
  *plhs = mxCreateSparse( m, n, nnz, mxREAL );
  A.val =       mxGetPr( *plhs );
  A.ind = (int*)mxGetIr( *plhs );
  A.col = (int*)mxGetJc( *plhs );
  return A;
}

/* =====================================================================
   Workspace.
   ====================================================================*/
typedef struct { /* BCLS-USER workspace. */
    matlab_Mat   Adense;
    matlab_spMat Asparse;
    mxArray      *argsRHS[4];
    mxArray      *Aprod_handle;
    mxArray      *Usolve_handle;
    double       *xTmp;         /* Used as temporary real workspace. */
} worksp;
#define FUNC_HANDLE argsRHS[0]
#define MODE        argsRHS[1]
#define IXvec       argsRHS[2]
#define Zvec        argsRHS[3]

/* ---------------------------------------------------------------------
   Inline max function.
   ---------------------------------------------------------------------*/
static int
imax(const int a, const int b) {
    return a > b ? a : b;
}

/* ---------------------------------------------------------------------
   Increase ctrl-c counter each time ctrl-c is pressed.
  ---------------------------------------------------------------------*/
static void
detect_ctrl_c_press(int sigID)
{
    ctrl_c_pressed += 1;
    if (ctrl_c_pressed == 1)
        mexPrintf( "     <ctrl-c> pressed. "
                   " BCLS will exit at end of major iteration.\n" );
    if (ctrl_c_pressed  > 1)
        mexPrintf( "     <ctrl-c> pressed (again!). "
                   " BCLS will exit at end of major iteration.\n" );
}

/* ---------------------------------------------------------------------
   CallBack.
   Return the number of times that ctrl-c was pressed.
   ---------------------------------------------------------------------*/
static int
CallBack( BCLS *ls, void *UsrWrk )
{
    return ctrl_c_pressed;
}

/* ---------------------------------------------------------------------
   pretty_printer
   Print to standard out.  Note that we can't rely on fprintf if
   Matlab is being run within its own GUI.
  ---------------------------------------------------------------------*/
static int
pretty_printer( void *io_file, char *msg)
{
    mexPrintf( msg );
    return 0;
}

/* ---------------------------------------------------------------------
   Usolve - calls the user-supplied preconditioning routine.
   ---------------------------------------------------------------------*/
static int
Usolve( int mode, int m, int n, int nix, int ix[],
        double v[], double w[], void *ptrA ) {

    worksp *Wrk = (worksp *)ptrA;
    int i, err;
    mxArray *plhs[1];
    double  *dTmp;
 
    mxAssert( mode == BCLS_PRECON_U    || 
              mode == BCLS_PRECON_Ut   ||
              mode == BCLS_PRECON_INIT ||
              mode == BCLS_PRECON_TERM,
              "Wrong mode." );

    /* Get the function handle to the user's Usolve routine. */
    Wrk->FUNC_HANDLE = Wrk->Usolve_handle;

    /* Set mode. */
    *(mxGetPr( Wrk->MODE )) = (double)mode;

    /* Set z(nix). */
    mxSetM( Wrk->Zvec, nix );
    if (mode == BCLS_PRECON_U)
        mxSetPr( Wrk->Zvec, w );
    else
        mxSetPr( Wrk->Zvec, v );

    /* Copy index information to Matlab's ix. */
    mxSetM( Wrk->IXvec, nix );
    dTmp = mxGetPr( Wrk->IXvec );
    for (i = 0; i < nix; i++) dTmp[i] = (double)(ix[i] + 1);

    /* Evaluate the user's aprod routine. */
    if (mode == BCLS_PRECON_INIT || mode == BCLS_PRECON_TERM) {
        signal(SIGINT, MATLAB_handler);
        err = mexCallMATLAB( 0, plhs, 4, Wrk->argsRHS, "feval" );
        signal(SIGINT, detect_ctrl_c_press);
        return err;
    }
    else {
        signal(SIGINT, MATLAB_handler);
        err = mexCallMATLAB( 1, plhs, 4, Wrk->argsRHS, "feval" );
        signal(SIGINT, detect_ctrl_c_press);
        if (err)
            return err;
    }

    /* Copy the result back into x or y and destroy LHS. */
    dTmp = mxGetPr( plhs[0] );
    if (mode == BCLS_PRECON_U)
        cblas_dcopy( nix, dTmp, 1, v, 1 );
    else if (mode == BCLS_PRECON_Ut)
        cblas_dcopy( nix, dTmp, 1, w, 1 );
    mxDestroyArray( plhs[0] );

    /* Successful exit. */
    return 0;
}

/* ---------------------------------------------------------------------
   AprodImplicit - matrix-vector products.

   If     mode == BCLS_PROD_A,  compute  y <- A *x, with  x  untouched;
   and if mode == BCLS_PROD_At, compute  x <- A'*y, with  y  untouched.
  ---------------------------------------------------------------------*/
static int
AprodImplicit( int mode, int m, int n, int nix, int ix[],
	       double x[], double y[], void *ptrA )
{
    worksp *Wrk = (worksp *)ptrA;
    int i, err;
    mxArray *plhs[1];
    double  *dTmp;
 
    mxAssert( mode == BCLS_PROD_A    || 
              mode == BCLS_PROD_At   ||
              mode == BCLS_PROD_INIT ||
              mode == BCLS_PROD_TERM,
              "Wrong mode." );

    /* Get the function handle to the user's Aprod routine. */
    Wrk->FUNC_HANDLE = Wrk->Aprod_handle;

    /* Set mode. */
    *(mxGetPr( Wrk->MODE )) = (double)mode;

    /* Set z(m or n) */
    if (mode == BCLS_PROD_A) {
        mxSetM ( Wrk->Zvec, n );
        mxSetPr( Wrk->Zvec, x );
    }
    else if (mode == BCLS_PROD_At) {
        mxSetM ( Wrk->Zvec, m );
        mxSetPr( Wrk->Zvec, y );
    }
    else { /* Empty matrix. */
        mxSetM ( Wrk->Zvec, 0 );
        mxSetPr( Wrk->Zvec, x );
    }
    
    /* Copy index information to Matlab's ix. */
    mxSetM( Wrk->IXvec, nix );
    dTmp = mxGetPr( Wrk->IXvec );
    for (i = 0; i < nix; i++) dTmp[i] = (double)(ix[i] + 1);

    /* Evaluate the user's aprod routine. */
    if (mode == BCLS_PROD_INIT || mode == BCLS_PROD_TERM) {
        signal(SIGINT, MATLAB_handler);
        err = mexCallMATLAB( 0, plhs, 4, Wrk->argsRHS, "feval" );
        signal(SIGINT, detect_ctrl_c_press);
        return err;
    }
    else {
        signal(SIGINT, MATLAB_handler);
        err = mexCallMATLAB( 1, plhs, 4, Wrk->argsRHS, "feval" );
        signal(SIGINT, detect_ctrl_c_press);
        if (err)
            return err;
    }

    /* Copy the result back into x or y and destroy LHS. */
    dTmp = mxGetPr( plhs[0] );
    if (mode == BCLS_PROD_A)
        cblas_dcopy( m, dTmp, 1, y, 1 );
    else
        cblas_dcopy( n, dTmp, 1, x, 1 );
    mxDestroyArray( plhs[0] );

    /* Successful exit. */
    return 0;
}

/* ----------------------------------------------------------------------
   AprodExplicitDense - dense matrix-vector products.

   If  mode == BCLS_PROD_A,  compute  y <- A *x, with  x  untouched;
   if  mode == BCLS_PROD_At, compute  x <- A'*y, with  y  untouched.
   ---------------------------------------------------------------------*/
static int
AprodExplicitDense( int mode, int m, int n, int nix, int ix[],
                    double x[], double y[], void *ptrA ) {

    worksp *Wrk     = (worksp *)ptrA;
    matlab_Mat *A   = &(Wrk->Adense);
    double *xTmp    = Wrk->xTmp;
    int i, j;

    mxAssert( mode == BCLS_PROD_A    || 
              mode == BCLS_PROD_At   ||
              mode == BCLS_PROD_INIT ||
              mode == BCLS_PROD_TERM,
              "Wrong mode." );
    mxAssert( A->m == m, "Wrong number of constraints." );
    mxAssert( A->n == n, "Wrong number of variables." );

    if (mode == BCLS_PROD_A) {
        
        /* Scatter x(ix) into xTmp. */
        for (j = 0; j < n; j++) xTmp[j] = 0.0;
        for (i = 0; i < nix; i++) {
            j = ix[i];
            xTmp[j] = x[j];
        }
        cblas_dgemv( CblasColMajor, CblasNoTrans,
                     m, n, 1.0, A->val, m, xTmp, 1, 0.0, y, 1 );
    }
    else if (mode == BCLS_PROD_At) {
        cblas_dgemv( CblasColMajor, CblasTrans,
                     m, n, 1.0, A->val, m, y, 1, 0.0, x, 1 );
    }    
    return 0;
}
/* ----------------------------------------------------------------------
   AprodExplicitSparse - sparse matrix-vector products.

   If  mode == BCLS_PROD_A,  compute  y <- A *x, with  x  untouched;
   if  mode == BCLS_PROD_At, compute  x <- A'*y, with  y  untouched.
   ---------------------------------------------------------------------*/
static int
AprodExplicitSparse( int mode, int m, int n, int nix, int ix[],
                     double x[], double y[], void *ptrA ) {

    worksp *Wrk = (worksp *)ptrA;
    matlab_spMat *A = &(Wrk->Asparse); 
    int    i, j, k, l;
    double aij, xj, sum;
    
    mxAssert( mode == BCLS_PROD_A    || 
              mode == BCLS_PROD_At   ||
              mode == BCLS_PROD_INIT ||
              mode == BCLS_PROD_TERM,
              "Wrong mode." );
    mxAssert( A->m == m, "Wrong number of constraints." );
    mxAssert( A->n == n, "Wrong number of variables." );
    
    if (mode == BCLS_PROD_A) {
	
	memset( y, 0, m * sizeof(double) );

	for (l = 0; l < nix; l++) {
	    j = ix[l];
	    xj = x[j];
	    if (xj == 0.0)
		; /* Relax. */
	    else
		for (k = A->col[j]; k < A->col[j+1]; k++) {
		    aij   = A->val[k];
		    i     = A->ind[k];
		    y[i] += aij * xj;
		}
	}
    }

    else if (mode == BCLS_PROD_At) {

	for (l = 0; l < nix; l++) {
	    j = ix[l];
	    sum = 0;
	    for (k = A->col[j]; k < A->col[j+1]; k++) {
		aij  = A->val[k];
		i    = A->ind[k];
		sum += aij * y[i];
	    }
	    x[j] = sum;
	}
    }

    /* Exit. */
    return 0;
}

/* ---------------------------------------------------------------------
   Main
   ---------------------------------------------------------------------*/
void mexFunction( int nlhs,       mxArray * plhs[],
		  int nrhs, const mxArray * prhs[] ) {

    enum {denseA, sparseA, operatorA} matrix_type;
    worksp       Wrk;
    int          i, j, err, buflen;
    double       *dTmp;
    void         *Adata;
    int (*Aprod)( int mode, int m, int n, int nix,
		  int ix[], double x[], double y[], void *UsrWrk );

    /* Install our own signal handler (and save Matlab's). */
    ctrl_c_pressed = 0;    
    MATLAB_handler = signal(SIGINT, detect_ctrl_c_press);

    /* Let the MEX interface handle any errors. */
    mexSetTrapFlag(1);

    /* -- Input 0-1: m, n. */
    int m = (int) assertScalar( prhs[0], "m" );
    int n = (int) assertScalar( prhs[1], "n" );

    /* Dimensions now known.  Initialize a BCLS problem. */
    BCLS *ls = bcls_create_prob( m, n );
    
    /* Input 3-5: b, bl, bu. */
    double *b   = assertMatrix( prhs[3], m, 1, "b"  );
    double *bl  = assertMatrix( prhs[4], n, 1, "bl" );
    double *bu  = assertMatrix( prhs[5], n, 1, "bu" );

    /* Input 6/Output 0: x. */
    plhs[0]     = mxDuplicateArray( prhs[6] );
    double *x   = assertMatrix( plhs[0], n, 1, "x"  );
    
    /* Input 7: c. */
    double *c = NULL;
    if (!mxIsEmpty(prhs[7])) c = assertMatrix( prhs[7], n, 1, "c" );

    /* Input 8: damp. */
    double damp = assertScalar( prhs[8], "damp" );

    /* Input 9: options. */
    const int numOpts = 8;
    mxArray *opt[numOpts];
    for (i = 0; i < numOpts; i++)
	opt[i] = mxGetFieldByNumber( prhs[9], 0, i );

    if (!mxIsEmpty(opt[0])) ls->print_level = (int)mxGetScalar(opt[0]);
    if (!mxIsEmpty(opt[1])) ls->proj_search = (int)mxGetScalar(opt[1]);
    if (!mxIsEmpty(opt[2])) ls->optTol      =      mxGetScalar(opt[2]);
    if (!mxIsEmpty(opt[3])) ls->itnMajLim   = (int)mxGetScalar(opt[3]);
    if (!mxIsEmpty(opt[4])) ls->itnMinLim   = (int)mxGetScalar(opt[4]);
    if (!mxIsEmpty(opt[5])) {
        bcls_set_usolve( ls, Usolve );
        Wrk.Usolve_handle = opt[5];
    }
    if (!mxIsEmpty(opt[6])) {
        /* Returns null if opt[6] isn't sensible. Fine for BCLS. */
        ls->minor_file = fopen(mxArrayToString(opt[6]), "w+");
    }
    if (!mxIsEmpty(opt[7])) ls->newton_step = (int)mxGetScalar(opt[7]);
    ls->CallBack = CallBack;

    /* -- Input 2: A. */
    if (mxIsNumeric(prhs[2])) {
        if (mxIsSparse(prhs[2])) matrix_type = sparseA;
        else                     matrix_type = denseA;
    }
    else
        matrix_type = operatorA;

    /* Allocate workspace for Implicit A and Usolve. */
    const int mnmax = imax( n, m );
    Wrk.MODE  = mxCreateDoubleScalar( 0            ); /* mode */
    Wrk.IXvec = mxCreateDoubleMatrix( n, 1, mxREAL ); /* ix   */
    Wrk.Zvec  = mxCreateDoubleMatrix( 1, 1, mxREAL ); /* z    */
    Wrk.xTmp  = mxCalloc( mnmax, sizeof(double) );    /* temp */
    
    /* The real pointer to Zvec will be overwritten in AprodImplicit.
       Save it and reset Zvec later to avoid memory leaks. */
    double *zPtr = mxGetPr(Wrk.Zvec);
    
    if (matrix_type == denseA) {
        Wrk.Adense.m     = m;
        Wrk.Adense.n     = n;
        Wrk.Adense.val   = assertMatrix( prhs[2], m, n, "A" );
        Aprod            = AprodExplicitDense;
    }
    else if (matrix_type == sparseA) {
	Wrk.Asparse.m    = m;
	Wrk.Asparse.n    = n;
	Wrk.Asparse.val  = assertSpMatrix( prhs[2], m, n, "A" );
	Wrk.Asparse.ind  = (int*) mxGetIr( prhs[2] );
	Wrk.Asparse.col  = (int*) mxGetJc( prhs[2] );
	Wrk.Asparse.nnz  = Wrk.Asparse.col[n];
	Aprod            = AprodExplicitSparse;
    }
    else {
	Aprod            = AprodImplicit;
        Wrk.Aprod_handle = (mxArray *)prhs[2];
    }
    
    bcls_set_print_hook( ls, stdout, pretty_printer );
    bcls_set_problem_data( ls,    /* The BCLS problem */
			   m,     /* Number of problem rows */
			   n,     /* Number of problem columns */
			   Aprod, /* The Mat-vec routine */
			   &Wrk,  /* Arbitrary data for the Mat-vec routine */
			   damp,  /* Damping parameter */
 			   x,     /* Solution vector */
			   b,     /* RHS vector */
			   c,     /* Linear term (may be NULL) */
			   bl,    /* Lower-bounds vector */
			   bu );  /* Upper-bounds vector */

    /* Call the main BCLS routine. */
    err = bcls_solve_prob( ls );
    
    /* Output 0: x. */
    /* - done - */
    
    /* Output 1: g. */
    plhs[1] = mxCreateDoubleMatrix( n, 1, mxREAL  );
    dTmp = mxGetPr( plhs[1] );
    for (i = 0; i < n; i++) dTmp[i] = ls->g[i];

    /* Output 2-5: dInf, jInf, majorIts, minorIts, status. */
    plhs[2] = mxCreateDoubleScalar( ls->soln_dInf );
    plhs[3] = mxCreateDoubleScalar( ls->soln_jInf );
    plhs[4] = mxCreateDoubleScalar( ls->itnMaj    );
    plhs[5] = mxCreateDoubleScalar( ls->itnMin    );
    plhs[6] = mxCreateDoubleScalar( ls->exit      );
    plhs[7] = mxCreateString( bcls_exit_msg( ls->exit ) );
    plhs[8] = mxCreateDoubleScalar( ls->stopwatch[BCLS_TIMER_TOTAL].total );
    
    /* Deallocate the BCLS problem. */
    err = bcls_free_prob( ls );
    mxAssert( err == 0 , "Couldn't deallocate BCLS context.");
    
    /* Restore the real pointer to Zvec, then deallocate all. */
    mxSetPr( Wrk.Zvec, zPtr );
    mxDestroyArray( Wrk.Zvec  );
    mxDestroyArray( Wrk.MODE  );
    mxDestroyArray( Wrk.IXvec );
    mxFree( Wrk.xTmp );

    /* Close files. */
    if (ls->minor_file) fclose( ls->minor_file );

    /* Reinstall Matlab's signal handler. */
    signal(SIGINT, MATLAB_handler);

    /* Exit. */
    return;
}
