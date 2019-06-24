/* =====================================================================
   $Revision: 290 $ $Date: 2007-03-04 22:01:54 -0800 (Sun, 04 Mar 2007) $
   Matlab MEX interface to iohb.c.

   25 Aug 05: Original version.
              Michael P. Friedlander
              mpf@cs.ubc.ca
              University of British Columbia
   ===================================================================== */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "mex.h"
#include "iohb.h"

/* =====================================================================
   Insert revision info into the object file.
   ===================================================================== */
static char *revisionid =
  "$Revision: 290 $"
  "$Date: 2007-03-04 22:01:54 -0800 (Sun, 04 Mar 2007) $";

/* =====================================================================
   Matlab helper routines.
   ===================================================================== */
typedef struct { /* Matlab full matrix. */
  int m, n;
  double *val;
} matlab_Mat;

typedef struct {  /* Compressed column (MATLAB) sparse format. */
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
  A.val = mxGetPr( *plhs );
  A.ind = (int*)mxGetIr( *plhs );
  A.col = (int*)mxGetJc( *plhs );
  return A;
}

/* ---------------------------------------------------------------------
   Main
   --------------------------------------------------------------------- */
void mexFunction( int nlhs,       mxArray * plhs[],
		  int nrhs, const mxArray * prhs[] ) {
   int    i, j;
   int    Nrhs     = 1;
   double *guess   = NULL;
   double *exact   = NULL;
   char   *Type    = "RRA";       /* (R)eal, (R)ectangular, (A)ssembled. */
   char   *RhsType = "F  ";       /* (F)ull, (G)uess, e(X)act soln. */
   char   Ptrfmt[] = "(13I6)";
   char   Indfmt[] = "(16I5)";
   char   Valfmt[] = "(3E26.18)";
   char   Rhsfmt[] = "(3E26.18)";

   double dTmp;
   char   msgbuf[60];

   FILE * out_file;
   char   hbf_file[100]; /* Name of Harwell-Boeing file. */
   char   opt_file[100]; /* Name of optional data file. */

   /* Gather the input data. */
   char *filename = assertString( prhs[0], "filename", 4 );
   char *title    = assertString( prhs[1], "title",    0 );
   char *key      = assertString( prhs[2], "key",      0 );

   matlab_spMat A;
   A.m    =        mxGetM ( prhs[3] );
   A.n    =        mxGetN ( prhs[3] );
   A.val  = assertSpMatrix( prhs[3], A.m, A.n, "A" );
   A.ind  = (int*) mxGetIr( prhs[3] );
   A.col  = (int*) mxGetJc( prhs[3] );
   A.nnz  = A.col[A.n];

   double *rhs = assertMatrix( prhs[4], A.m, 1, "b"  );
   double *bl  = assertMatrix( prhs[5], A.n, 1, "bl" );
   double *bu  = assertMatrix( prhs[6], A.n, 1, "bu" );
   double *c   = assertMatrix( prhs[7], A.n, 1, "c"  );
   double *x   = assertMatrix( prhs[8], A.n, 1, "x"  );
   double *odata[4] = {bl, bu, c, x};
   
   /* Write out A and rhs to the Harwell-Boeing file. */
   sprintf( hbf_file, "%s.%s", filename, "hbf" );
   
   writeHB_mat_double( hbf_file, A.m, A.n, A.nnz, A.col, A.ind, A.val,
		       Nrhs, rhs, guess, exact,
		       title, key, Type,
		       Ptrfmt, Indfmt, Valfmt, Rhsfmt, "F" );

   /* Write out bl and bu to the bounds file. */
   sprintf( opt_file, "%s.%s", filename, "opt" );
   out_file = fopen( opt_file, "w" );

   for (j = 0; j < 4; j++ ) {
       for (i = 0; i < A.n; i++) {
	   dTmp = odata[j][i];
	   if ( isnormal( dTmp ) ) 
	       ; /* Success.  Relax. */
	   if ( !isfinite( dTmp ) ) {
	       /* Value is infinite. */
	       if ( dTmp < 0 ) dTmp = -1e20;
	       else            dTmp =  1e20;
	   }
	   else if ( isnan( dTmp ) ) {
	       /* A NaN was detected.  Exit. */
	       sprintf( msgbuf, "NaN detected in position %d.\n", i+1);
	       goto error_exit;
	   }
	   fprintf( out_file, "%23.16e\n", dTmp );
       }
   }
   fclose( out_file );
   
   /* Clean exit. */
   return;
   
   /* Error exit. */
 error_exit:
   fclose( out_file );
   mexErrMsgTxt( msgbuf );

}
