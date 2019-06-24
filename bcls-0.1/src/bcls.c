/* bcls.c
   $Revision: 282 $ $Date: 2006-12-17 17:38:00 -0800 (Sun, 17 Dec 2006) $

   ---------------------------------------------------------------------
   This file is part of BCLS (Bound-Constrained Least Squares).

   Copyright (C) 2006 Michael P. Friedlander, Department of Computer
   Science, University of British Columbia, Canada. All rights
   reserved. E-mail: <mpf@cs.ubc.ca>.

   BCLS is free software; you can redistribute it and/or modify it
   under the terms of the GNU Lesser General Public License as
   published by the Free Software Foundation; either version 2.1 of the
   License, or (at your option) any later version.

   BCLS is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General
   Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with BCLS; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
   USA
   ---------------------------------------------------------------------
*/
/*!
  \file
  BCLS user-callable library routines.
*/

#include <cblas.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <assert.h>
#include <setjmp.h>

#include "bcls.h"
#include "bclib.h"
#include "bcsolver.h"
#include "bcversion.h"

/*!

  \brief Malloc wrapper.  Private to this file.

  Pointer to the newly allocated memory,
  or NULL if malloc returned an eror.
  
  \param[in]  len   Number of elements needed
  \param[in]  size  Memory needed for each element.
  \param[in]  who   Short character description of the memory needed.
  
  \return Pointer to the newly allocated memory.

*/
static void *
xmalloc (int len, size_t size, char * who)
{
    register void *value = malloc(len * size);
    if (value == NULL)
	fprintf( stderr, "Not enough memory to allocate %s.\n", who );
    return value;
}

/*!

  \brief Create a BCLS problem instance.

  This routine create a BCLS problem.  It *must* be called before any
  other BCLS routine.  It will create and initialize a new BCLS
  problem with enough space to accomodate the specified problem.  The
  problem is then initialized using bcls_init_prob.

  In order to re-initialize the BCLS problem instance, but not
  deallocate memory that's already sufficient for an equal or smaller
  size problem, use bcls_init_prob.

  \param[in] mmax  Maximum number of rows    in any A.
  \param[in] nmax  Maximum number of columns in any A.

  \return Return a pointer to a new BCLS problem.  If the return is
  NULL, then the initialization failed.

*/
BCLS *
bcls_create_prob( int mmax, int nmax )
{
    int mnmax = imax( nmax, mmax );

    assert( mmax > 0 && nmax > 0 );
    
    // Allocate the problem.
    BCLS *ls = (BCLS *)xmalloc(1, sizeof(BCLS), "ls" );

    // Check if the problem has been successfully allocated.
    if (ls == NULL) {
	fprintf( stderr, "XXX Could not allocate a BCLS problem.\n");
	return ls;
    }

    // -----------------------------------------------------------------
    // Allocate workspace vectors large enough to accomdate the problem.
    // -----------------------------------------------------------------

    // Record the maximum alloted workspace.
    ls->mmax = mmax;
    ls->nmax = nmax;

    // Residual: r(m+n).
    ls->r = (double *)xmalloc( mmax+nmax, sizeof(double), "r" );
    if (ls->r == NULL) goto error;

    // Gradient: g(n).
    ls->g = (double *)xmalloc( nmax, sizeof(double), "g" );
    if (ls->g == NULL) goto error;

    // Search direction, full space: dx(n).
    ls->dx = (double *)xmalloc( nmax, sizeof(double), "dx" );
    if (ls->dx == NULL) goto error;

    // Search direction, subspace: dxFree(n).
    ls->dxFree = (double *)xmalloc( nmax, sizeof(double), "dxFree" );
    if (ls->dxFree == NULL) goto error;

    // Step to each breakpoint: aBreak(n).
    ls->aBreak = (double *)xmalloc( nmax, sizeof(double), "aBreak" );
    if (ls->aBreak == NULL) goto error;

    // Indices of each breakpoint: iBreak(n).
    ls->iBreak = (int *)xmalloc( nmax, sizeof(int), "iBreak" );
    if (ls->iBreak == NULL) goto error;

    // Variable indices: ix(n).
    ls->ix = (int *)xmalloc( nmax, sizeof(int), "ix" );
    if (ls->ix == NULL) goto error;

    // Workspace: wrk_v( max(n, m) ).
    ls->wrk_u = (double *)xmalloc( mnmax, sizeof(double), "wrk_u" );
    if (ls->wrk_u == NULL) goto error;

    // Workspace: wrk_v( max(n, m) ).
    ls->wrk_v = (double *)xmalloc( mnmax, sizeof(double), "wrk_v" );
    if (ls->wrk_v == NULL) goto error;

    // Workspace: wrk_w( max(n, m) ).
    ls->wrk_w = (double *)xmalloc( mnmax, sizeof(double), "wrk_w" );
    if (ls->wrk_w == NULL) goto error;

    // -----------------------------------------------------------------
    // Initialize this new problem instance.
    // -----------------------------------------------------------------
    bcls_init_prob( ls );

    // -----------------------------------------------------------------
    // Exits.
    // -----------------------------------------------------------------

    // Successfull exit.
    return ls;

 error:
    // Unsuccessful exit.
    bcls_free_prob( ls );
    return ls;
}

/*!

  \brief Initialize a BCLS problem.
  
  Initialize a BCLS problem.  You can call this routine to "reset" a
  BCLS problem, i.e., reset all parameter values.  Note that it is
  automatically called by bcls_create_prob.

  \param[in,out] ls  BCLS problem context.

*/
void
bcls_init_prob( BCLS *ls )
{
    // Some quick error checking.
    assert( ls->mmax > 0 && ls->nmax > 0 );

    // Check if the problem has been successfully allocated.
    if (ls == NULL) {
	fprintf( stderr, "XXX The BCLS problem is NULL.\n");
	return;
    }
    
    // Initialize the problem structure.
    ls->print_info    =  NULL;
    ls->print_hook    =  NULL;
    ls->fault_info    =  NULL;
    ls->fault_hook    =  NULL;
    ls->Aprod         =  NULL;
    ls->Usolve        =  NULL;
    ls->CallBack      =  NULL;
    ls->UsrWrk        =  NULL;
    ls->anorm         =  NULL;
    ls->print_level   =  1;
    ls->proj_search   =  BCLS_PROJ_SEARCH_EXACT;
    ls->newton_step   =  BCLS_NEWTON_STEP_LSQR;
    ls->minor_file    =  NULL;
    ls->itnMaj        =  0;
    ls->itnMajLim     =  5  * (ls->nmax);
    ls->itnMin        =  0;
    ls->itnMinLim     =  10 * (ls->nmax);
    ls->nAprodT       =  0;
    ls->nAprodF       =  0;
    ls->nAprod1       =  0;
    ls->nUsolve       =  0;
    ls->m             =  0;
    ls->n             =  0;
    ls->unconstrained =  0;
    ls->damp          =  0.0;
    ls->damp_min      =  1.0e-4;
    ls->exit          =  BCLS_EXIT_UNDEF;
    ls->soln_rNorm    = -1.0;
    ls->soln_dInf     =  0.0;
    ls->soln_jInf     = -1;
    ls->soln_stat     =  BCLS_SOLN_UNDEF;
    ls->optTol        =  1.0e-6;
    ls->conlim        =  1.0 / ( 10.0 * sqrt( DBL_EPSILON ) );
    ls->mu            =  1.0e-2;
    ls->backtrack     =  1.0e-1;
    ls->backtrack_limit = 10;

    // Initialize timers.
    bcls_timer( &(ls->stopwatch[BCLS_TIMER_TOTAL] ),  BCLS_TIMER_INIT );
    ls->stopwatch[BCLS_TIMER_TOTAL].name = "Total time";

    bcls_timer( &(ls->stopwatch[BCLS_TIMER_APROD] ),  BCLS_TIMER_INIT );
    ls->stopwatch[BCLS_TIMER_APROD].name = "Total time for Aprod";

    bcls_timer( &(ls->stopwatch[BCLS_TIMER_USOLVE] ),  BCLS_TIMER_INIT );
    ls->stopwatch[BCLS_TIMER_USOLVE].name = "Total time for Usolve";

    bcls_timer( &(ls->stopwatch[BCLS_TIMER_LSQR] ),  BCLS_TIMER_INIT );
    ls->stopwatch[BCLS_TIMER_LSQR].name = "Total time for LSQR";

    // Initialize constants.
    ls->eps         = DBL_EPSILON;
    ls->eps2        = pow( DBL_EPSILON, 0.50 );
    ls->eps3        = pow( DBL_EPSILON, 0.75 );
    ls->epsx        = ls->eps3;
    ls->epsfixed    = DBL_EPSILON;
    ls->BigNum      = BCLS_INFINITY;

    return;
}

/*!

  \brief Free and reinitialize all resources in an existing BCLS problem.

  \return
  - 0: no errors
  - 1: the BCLS problem is NULL (ie, does not exist).

*/
int
bcls_free_prob( BCLS *ls )
{
    // Report an error if it's already NULL.
    if (ls == NULL) return 1;

    // Deallocate the workspace.
    assert( ls->r       != NULL ); free( ls->r       );
    assert( ls->g       != NULL ); free( ls->g       );
    assert( ls->dx      != NULL ); free( ls->dx      );
    assert( ls->dxFree  != NULL ); free( ls->dxFree  );
    assert( ls->aBreak  != NULL ); free( ls->aBreak  );
    assert( ls->iBreak  != NULL ); free( ls->iBreak  );
    assert( ls->ix      != NULL ); free( ls->ix      );
    assert( ls->wrk_u   != NULL ); free( ls->wrk_u   );
    assert( ls->wrk_v   != NULL ); free( ls->wrk_v   );
    assert( ls->wrk_w   != NULL ); free( ls->wrk_w   );
    
    // Deallocate the BCLS problem.
    free( ls );
    
    return 0;
}

/*!

  
  \brief Install a print-hook routine.

  This routine installs a user-defined print-hook routine.
  
  The parameter info is a transit pointer passed to the hook routine.
  
  The parameter hook is an entry point to the user-defined print-hook
  routine. This routine is called by the routine "print" every time an
  informative message should be output. The routine "print" passes to
  the hook routine the transit pointer info and the character string
  msg, which contains the message. If the hook routine returns zero,
  the routine print prints the message in an usual way. Otherwise, if
  the hook routine returns non-zero, the message is not printed.
  
  In order to uninstall the hook routine the parameter hook should be
  specified as NULL (in this case the parameter info is ignored).

  \param[in,out] ls    BCLS problem context.
  \param[in]     info  Transit pointer passed to the hook routine.
  \param[in]     hook  Pointer to the print-hook routine.

*/
void
bcls_set_print_hook( BCLS *ls,
		     void *info,
		     int (*hook)(void *info, char *msg))
{
    ls->print_info = info;
    ls->print_hook = hook;
    return;
}

/*!

  \brief Install a fault-hook routine.

  This routine installs a user-defined fault-hook routine.

  The parameter info is a transit pointer passed to the hook routine.

  The parameter "hook" is an entry point to the user-defined
  fault-hook routine. This routine is called by the routine
  "bcls_fault" every time an error message should be output. The
  routine "bcls_fault" passes to the hook routine the transit pointer
  info and the character string msg, which contains the message. If
  the hook routine returns zero, the routine print prints the message
  in an usual way. Otherwise, if the hook routine returns non-zero,
  the message is not printed.

  In order to uninstall the hook routine the parameter hook should be
  specified as NULL (in this case the parameter info is ignored).

  \param[in,out] ls    BCLS problem context.
  \param[in]     info  Transit pointer passed to the hook routine.
  \param[in]     hook  Pointer to the fault-hook routine.

*/
void
bcls_set_fault_hook( BCLS *ls,
		     void *info,
		     int (*hook)(void *info, char *msg))
{
    ls->fault_info = info;
    ls->fault_hook = hook;
    return;
}

/*!

  \brief Give BCLS access to set of column weights.

  Define a set of column weights.  BCLS will use these to scale the
  steepest-descent steps.

  \param[in,out] ls     BCLS problem context.
  \param[in]     anorm  Array of columns norms of A.

*/
void
bcls_set_anorm( BCLS *ls,
                double anorm[] )
{
    ls->anorm = anorm;
    return;
}

/*!

  \brief Install a user-defined preconditioning routine.
 
  Replace each subproblem
  \f[
  \def\minimize{\displaystyle\mathop{\hbox{minimize}}}
  \minimize_x \| Ax - b \|
  \f]
  with
  \f[
  \def\minimize{\displaystyle\mathop{\hbox{minimize}}}
  \minimize_y \| AU^{-1}y - b \|
  \quad\mbox{with}\quad
  Ux = y.
  \f]

  \param[in,out]  ls     BCLS problem context.
  \param[in]      Usolve Pointer to the user's preconditioning routine.

*/
void
bcls_set_usolve( BCLS *ls,
                 int (*Usolve)( int mode, int m, int n, int nix,
                                int ix[], double v[], double w[],
                                void *UsrWrk ) )
{
    assert( Usolve != NULL );
    ls->Usolve = Usolve;
    return;
}


/*!

  \brief Compute the column norms of A.

  Compute the columns norms of A.  Each column of A is generated by
  multiplying A by a unit vector.  The two-norm squared of each
  column is stored in aprod.
  
  This routine must be called *after* bcls_create_prob.  (As with any
  other BCLS routine.)

  \param[in]      ls      BCLS problem context.
  \param[in]      m       Number of rows in A.
  \param[in]      n       Number of columns in A.
  \param[in]      Aprod   User's matrix-product routine.
  \param[in,out]  UsrWrk  Pointer to user's workspace (could be NULL).
  \param[out]     anorm   Vector of column norms of A.
  
  \return Returns any error codes returned by the user's aprod
  routine.

*/
int
bcls_compute_anorm( BCLS *ls, int n, int m,
                    int (*Aprod)
                    ( int mode, int m, int n, int nix,
                      int ix[], double x[], double y[], void *UsrWrk ),
                    void *UsrWrk,
                    double anorm[] )
{    
    // Make sure that the problem instance has been created.
    assert( ls    != NULL );
    assert( anorm != NULL );

    int j;
    int err;
    const int mode = 1;         // Always:  aj = A*ej.
    const int nix  = 1;         // Only one column is ever needed.
    int    *ix = ls->ix;
    double *e  = ls->wrk_u;
    double *aj = ls->wrk_v;

    bcls_dload( n, 1.0, e, 1 );
    
    for (j = 0; j < n; j++) {

        // Multiply A times the jth unit vector.
        ix[0] = j;
        
        // aj <- A * ej.
        err = Aprod( mode, m, n, nix, ix, e, aj, UsrWrk );

        // Exit if Aprod returned an error.
        if (err) break;

        // Compute the norm of aj and store it in anorm.
        anorm[j] = cblas_dnrm2( m, aj, 1 );

        // Make sure that the column norm is too small.
        anorm[j] = fmax( BCLS_MIN_COLUMN_NORM, anorm[j] );
    }

    return err;
}

/*!

  \brief Load the problem into the BCLS data structure.

  The routiens loads the complete problem description into a BCLS
  problem instance.  No work is really being done: only the various
  structure elements are being set to point to the right places.

  \param[in,out]  ls      BCLS problem context.
  \param[in]      m       Number of rows in A.  Note: m <= mmax.
  \param[in]      n       Number of columns in A.  Note: n <= nmax.
  \param[in]      Aprod   Hook to user's matrix-vector product routine.
  \param[in,out]  UsrWrk  Transit pointer passed directly to Aprod/Usolve.
  \param[in]      damp    Regularization parameter.
  \param[in,out]  x       The starting point.
  \param[in]      b       The RHS vector.
  \param[in]      c       Defines a linear term (set to NULL if one doesn't exist).
  \param[in]      bl      Lower bounds on x.
  \param[in]      bu      Upper bounds on x.

*/
void
bcls_set_problem_data( BCLS *ls, int m, int n,
		       int (*Aprod)( int mode, int m, int n, int nix,
				     int ix[], double x[], double y[],
				     void *UsrWrk ),
		       void *UsrWrk, double damp,
		       double x[], double b[], double c[],
		       double bl[], double bu[] )
{
    assert( m <= ls->mmax );
    assert( n <= ls->nmax );

    assert( Aprod  != NULL ); ls->Aprod  = Aprod;
    assert( m      >= 0    ); ls->m      = m;
    assert( n      >= 0    ); ls->n      = n;
                              ls->UsrWrk = UsrWrk;
    assert( damp   >= 0.0  ); ls->damp   = damp;
    assert( x      != NULL ); ls->x      = x;
    assert( b      != NULL ); ls->b      = b;
                              ls->c      = c;
    assert( bl     != NULL ); ls->bl     = bl;
    assert( bu     != NULL ); ls->bu     = bu;

    return;
}

/*!

  \brief Return an exit message.

  \param[in] flag  The code of the error.

  \return Error messsage.

*/
char *
bcls_exit_msg( int flag )
{
    char *msg;

    if      (flag==BCLS_EXIT_CNVGD)  msg="Optimal solution found";
    else if (flag==BCLS_EXIT_MAJOR)  msg="Too many major iterations";
    else if (flag==BCLS_EXIT_MINOR)  msg="Too many minor iterations";
    else if (flag==BCLS_EXIT_UNBND)  msg="Found direction of infinite descent";
    else if (flag==BCLS_EXIT_INFEA)  msg="Bounds are inconsistent";
    else if (flag==BCLS_EXIT_APROD)  msg="Aprod requested immediate exit";
    else if (flag==BCLS_EXIT_USOLVE) msg="Usolve requested immediate exit";
    else if (flag==BCLS_EXIT_CALLBK) msg="CallBack requested immediate exit";
    else                             msg="Undefined exit";

    return msg;
}

/*!

  \brief Solve the current problem instance.

  \return
   -  #BCLS_EXIT_CNVGD (0) - Successful exit.  Found an optimal solution.
   -  #BCLS_EXIT_MAJOR     - Too many major iterations.
   -  #BCLS_EXIT_MINOR     - Too many inner iterations.
   -  #BCLS_EXIT_UNBND     - Found direction of infinite descent.
   -  #BCLS_EXIT_INFEA     - Bounds are inconsistent.
   -  #BCLS_EXIT_APROD     - Aprod  requested immediate exit.
   -  #BCLS_EXIT_USOLVE    - Usolve requested immediate exit.

*/
int
bcls_solve_prob( BCLS *ls )
{
    int err;
    int jpInf;
    const int preconditioning = ls->Usolve != NULL;
    double pInf, timeTot;
    double bNorm;  // Norm of the RHS (or 1.0, which ever is bigger).
    BCLS_timer watch;

    // Print the banner.
    //PRINT1(" ----------------------------------------------------\n");
    //PRINT1(" BCLS -- Bound-Constrained Least Squares, Version %s\n",
    //      bcls_version_info() );
    //PRINT1(" Compiled on %s\n",  bcls_compilation_info() );
    //PRINT1(" ----------------------------------------------------\n");

    // -----------------------------------------------------------------
    // Print parameters and diagnostic information.
    // -----------------------------------------------------------------
    bcls_print_params( ls );

    // Examine the bounds and print a summary about them.
    // Also check if the problem is feasible!
    err = bcls_examine_bnds( ls, ls->n, ls->bl, ls->bu );
    if (err) {
	ls->exit = err;
        goto direct_exit;
    }

    // Set the return for longjmp.  First call to setjmp returns 0.
    err = setjmp(ls->jmp_env);
    if (err) {
        ls->exit = err;
        goto direct_exit;
    }

    // Compute and print some statistics about the column scales.
    bcls_examine_column_scales( ls, ls->anorm );

    // Print a log header.
    //    PRINT1("\n %5s  %5s    %9s  %9s  %11s  %9s %7s %7s %7s\n",
    //	   "Major","Minor","Residual","xNorm","Optimal",
    //	   ls->newton_step == BCLS_NEWTON_STEP_LSQR
    //	   ? "LSQR" : "CGLS",
    //	   "nSteps", "Free", "OptFac");

    // -----------------------------------------------------------------
    // Call the BCLS solver.
    // -----------------------------------------------------------------
    bcls_timer( &(ls->stopwatch[BCLS_TIMER_TOTAL]), BCLS_TIMER_START );

    bcls_solver( ls, ls->m, ls->n, &bNorm,
		 ls->x, ls->b, ls->c, ls->bl, ls->bu, ls->r, ls->g,
		 ls->dx, ls->dxFree, ls->ix,
		 ls->aBreak, ls->iBreak, ls->anorm );

    bcls_timer( &(ls->stopwatch[BCLS_TIMER_TOTAL]), BCLS_TIMER_STOP );


    // =================================================================
    // Exits.
    // =================================================================
 direct_exit:
    // -----------------------------------------------------------------
    // Check feasibility and print a solution summary (iterations, etc.)
    // -----------------------------------------------------------------
    bcls_primal_inf( ls->n, ls->x, ls->bl, ls->bu, &pInf, &jpInf );

    //PRINT1("\n");
    //    PRINT1(" BCLS Exit %4d -- %s\n",ls->exit, bcls_exit_msg(ls->exit));
    //PRINT1("\n");
    //PRINT1(" No. of iterations      %8d %7s", ls->itnMin, "");
    //PRINT1(" Objective value       %17.9e\n", ls->soln_obj  );
    //PRINT1(" No. of major iterations%8d %7s", ls->itnMaj, "");
    //PRINT1(" Optimality           (%7d) %8.1e\n",
    //	   ls->soln_jInf,ls->soln_dInf);
    //PRINT1(" No. of calls to Aprod  %8d", ls->nAprodT);

    //if (ls->proj_search == BCLS_PROJ_SEARCH_EXACT)
    //   PRINT1(" (%5d)", ls->nAprod1 );
    //else
    //  PRINT1("  %5s ", "");
    // PRINT1(" Norm of RHS            %16.1e\n", bNorm);
    // if (pInf > ls->eps)
    //  PRINT1(" Feasibility          (%7d) %7.1e\n", jpInf, pInf);
    //if (preconditioning)
    //  PRINT1(" No. of calls to Usolve %8d\n", ls->nUsolve);
    //PRINT1("\n");

    // -----------------------------------------------------------------
    // Print timing statistics.
    // -----------------------------------------------------------------
    //timeTot = fmax(ls->eps,ls->stopwatch[BCLS_TIMER_TOTAL].total);

    //watch   = ls->stopwatch[BCLS_TIMER_TOTAL];
    //timeTot = fmax(1e-6, watch.total); // Safeguard against timeTot = 0.

    //PRINT1(" %-25s %5.1f (%4.2f) secs\n",
    //	   watch.name, watch.total, 1.0 );

    //watch = ls->stopwatch[BCLS_TIMER_LSQR];
    //PRINT1(" %-25s %5.1f (%4.2f) secs\n",
    //	   watch.name, watch.total, watch.total / timeTot );

    //    watch = ls->stopwatch[BCLS_TIMER_APROD];
    //PRINT1(" %-25s %5.1f (%4.2f) secs\n",
    //	   watch.name, watch.total, watch.total / timeTot );

    // Only print time for Usolve is user provided this routine.
    //if (preconditioning) {
      //  watch = ls->stopwatch[BCLS_TIMER_USOLVE];
      //PRINT1(" %-25s %5.1f (%4.2f) secs\n",
      //       watch.name, watch.total, watch.total / timeTot );
    //}

    return (ls->exit + err);
}
