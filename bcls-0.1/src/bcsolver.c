/* bcsolver.c
   $Revision: 283 $ $Date: 2006-12-17 20:27:56 -0800 (Sun, 17 Dec 2006) $

   ----------------------------------------------------------------------
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
   ----------------------------------------------------------------------
*/
/*!
   \file
   The main algorithm for the BCLS solver.
*/

#include <cblas.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "bcls.h"
#include "bclib.h"
#include "bccgls.h"
#include "bclsqr.h"
#include "bcsolver.h"

// =====================================================================
// Private (static) function declarations.
// =====================================================================

// Compute a Newton step on the free variables, dx(free).
static int
bcls_newton_step( BCLS *ls, int m, int nFree, int ix[], double damp,
		  int itnLim, double tol, double dxFree[], double x[],
		  double c[], double r[], int *itns, double *rTol );

// Find the minimizer along the projected search direction.
static int
bcls_proj_search( BCLS *ls, int m, int n, double x[], double dx[],
		  double f, double g[], double aBreak[], int iBreak[],
		  int ix[], double bl[], double bu[],
		  double ex[], double Ad[], double Ae[], int *totHits );

// Simple backtracking to satisfy an Armijo condition.
static int
bcls_proj_backtrack( BCLS *ls, int m, int n, double x[], double dx[],
		     double f, double g[], double bl[], double bu[],
		     double s[], int ix[], double As[], int *nSteps );


/*!

  \brief Implementation of the BCLS projected Newton/gradient search.

  \param[in,out]  ls     The BCLS problem context.
  \param[in]      m      Number of rows in A.
  \param[in]      n      Number of columns in A.
  \param[in,out]  bNorm  Norm of the RHS.
  \param[in,out]  x      The current iterate.
  \param[in]      b      The RHS vector.
  \param[in]      c      The linear-term vector.
  \param[in]      bl     Lower bounds on x.
  \param[in]      bu     Upper bounds on x.
  \param[in,out]  r      The residual vector r = b - Ax.
  \param[in,out]  g      The objective gradient.
  \param[in,out]  dx     The step direction in the full space.
  \param[in,out]  dxFree The step direction on the free variables.
  \param[in,out]  ix     Indices of the free variables.
  \param[in,out]  aBreak A list of breakpoints for a search direction.
  \param[in,out]  iBreak The indices of the breakpoints.
  \param[in,out]  anorm  Column norms of A.

*/
void
bcls_solver( BCLS *ls, int m, int n, double *bNorm,
	     double x[], double b[], double c[], double bl[], double bu[],
	     double r[], double g[], double dx[], double dxFree[],
	     int ix[], double aBreak[], int iBreak[], double anorm[] )
{
    int
        err = 0,       // Error indicator.
        j,             // Misc. counter.
        jInf,          // Variable with largest dual infeasibility.
        itnLimk,       // No. of allowed itns for kth subproblem.
        itnsk = 0,     // No. of iterations performed on kth subproblem.
        nFree,         // No. of variables floating between their bounds.
        nSteps = 0;    // No. of projected search steps.

    double
        optk = 0.0,    // kth     subproblem achieved optimality.
        optkm1 = 0.0,  // (k-1)th ...
        subTol,        // Required subproblem optimality.
        rTol = 0.1,    // Inexact Newton optimality reduction.
        dInf,          // Dual infeasibility.
        rNorm,         // Norm of the residual.
        xNorm,         // Norm of the current solution estimate.
        f;             // Objective value = 1/2 * rNorm^2.

    const int scaled_steepest = anorm != NULL;
    const int linear     = c != NULL;
    const int damped     = ls->damp > 0.0;
    const double damp    = ls->damp;
    const double damp2   = damp*damp;
    const double rTolMin = ls->optTol;

    // Initialize various variables.
    ls->itnMaj = 0;    // Reset the major iteration counter.
    ls->itnMin = 0;    // Reset the minor iteration counter.

    // Push x into the box.
    bcls_mid( n, x, bl, bu );

    // Convergence will be normalized by ||A'b|| + ||c||.
    for (j = 0; j < n; j++) ix[j] = j;
    bcls_aprod( ls, BCLS_PROD_At, n, ix, g, b );
    const double AbNorm = cblas_dnrm2( n, g, 1 );             // ||A'b||
    const double cnorm  = (c) ? cblas_dnrm2( n, c, 1 ) : 0.0; // ||  c||
    *bNorm = cblas_dnrm2( m, b, 1 );                          // ||  b||

    // -----------------------------------------------------------------
    // Start BCLS major iterations.
    // -----------------------------------------------------------------
    while (1) {

        // Need to include contribution for *all* cols of A.
        for (j = 0; j < n; j++) ix[j] = j;

        // Evaluate the (unregularized) residual:
        // r(1:m) = b - Ax
        bcls_aprod( ls, BCLS_PROD_A, n, ix, x, r );  // r(1:m) =  Ax
        cblas_daxpy( m, -1.0, b, 1, r, 1 );          // r(1:m) = -b + Ax
        cblas_dscal( m, -1.0, r, 1 );                // r(1:m) =  b - Ax

        // Evaluate the gradient:
        // g = A'(Ax - b)  +  damp^2 x = -A'r + damp^2 x.
        bcls_aprod( ls, BCLS_PROD_At, n, ix, g, r ); // g =  A'r
        cblas_dscal( n, -1.0, g, 1 );                // g = -A'r
        if (damped)
            cblas_daxpy( n, damp2, x, 1, g, 1 );     // g = g + damp^2 x
        if (linear)
            cblas_daxpy( n, 1.0, c, 1, g, 1 );       // g = g + c
        
        // Norm of the current solution and residual.
        xNorm = cblas_dnrm2( n, x, 1 );              //  || x ||
        rNorm = cblas_dnrm2( m, r, 1 );              //  || r ||
            
        // Objective value:  f = 1/2 ||r||^2 + 1/2 damp^2 ||x||^2.
        f = 0.5 * rNorm * rNorm;
        if (damped)
            f = f + 0.5 * damp2 * xNorm * xNorm;
        if (linear)
            f = f + cblas_ddot( n, c, 1, x, 1 );

        // -------------------------------------------------------------
        // Evaluate optimality, active set, and exit conditions.
        // -------------------------------------------------------------

        // Store free variable indices in ix[:nFree].
        nFree = bcls_free_vars( ls, n, ix, x, g, bl, bu );

        // Evaluate optimality of the current iterate x.
        bcls_dual_inf( n, x, g, bl, bu, &dInf, &jInf );

        // Test for exit conditions.  Acted on below.
        if ( dInf <= ls->optTol * (1.0 + AbNorm + cnorm) ) {
            ls->exit      = BCLS_EXIT_CNVGD;
            ls->soln_stat = BCLS_SOLN_OPTIM;    // Optimal soln found.
        }
        else if ( ls->itnMaj >= ls->itnMajLim ) // Too many majors.
            ls->exit      = BCLS_EXIT_MAJOR;
        else if ( ls->itnMin >= ls->itnMinLim ) // Too many minors.
            ls->exit      = BCLS_EXIT_MINOR;
        else if ( bcls_callback( ls ) )         // Exit requested by user.
            ls->exit      = BCLS_EXIT_CALLBK;

        // -------------------------------------------------------------
        // Print out the "kth" iteration log.
        // -------------------------------------------------------------
        /*PRINT1(" %5d  %5d  %11.4e  %9.2e  %11.4e  %9.2e %7d %7d",
               ls->itnMaj, itnsk, rNorm, xNorm, dInf, optk, nSteps, nFree );
        
        // Check if rTol was "tight enough" for the last major.
        if ( ls->itnMaj <  1 ) {
            PRINT1(" %7.0e", rTol);
        }
        else if ( optk  <= .2 * optkm1 )
            ; // Relax. Optimality has made a reasonable reduction.
        else if ( rTol  <= rTolMin ) {
              // Inexact Newton factor is alrady small.
            PRINT1(" %7.0e", rTol);
        }
        else if ( itnsk <= 1 ) {
            rTol /= 10.0;
            PRINT1(" %7.0e", rTol);
        }

        // Finish the log (eg, print a new-line).
        PRINT1("\n");
	*/
        // Execute the exit condition (possibly set above).
        if ( ls->exit != BCLS_EXIT_UNDEF ) break;

        // =============================================================
        // ----  New major iteration starts here ----
        // =============================================================
        ls->itnMaj++;

        // Zero out the step direction: dx = 0.
        bcls_dload( n, 0.0, dx, 1 );

        // -------------------------------------------------------------
        // Compute a Newton step on the free variables: dxFree.
        // -------------------------------------------------------------
        if (nFree) {

            // Find out the number of remaining allowed minors.
            itnLimk = ls->itnMinLim - ls->itnMin;

            // Save the previous subproblem optimality value and
            // compute the required optimality for this next
            // subproblem.
            optkm1 = optk;
            if (ls->unconstrained)
                subTol = rTol * ls->optTol;
            else
                subTol = rTol;// * bcls_vec_dnorm2( nFree, ix, g );
            
            // Note that ls->dx is used as workspace.
            err = bcls_newton_step( ls, m, nFree, ix, damp, itnLimk,
                                    subTol, dxFree, x, c, r,
                                    &itnsk, &optk );
            
            // Accumulate the minor iteration count.
            ls->itnMin += itnsk;
        }
        
        // -------------------------------------------------------------
        // Load steepest descent into the fixed variables: dx(fixed).
        // -------------------------------------------------------------
        // Scale steepest descent if user provided column scales.
        // Note that this needs to be done *after* bcls_newton_step.
        cblas_dcopy( n, g, 1, dx, 1 );      //  dx <-  g
        cblas_dscal( n, -1.0, dx, 1 );      //  dx <- -g

        if ( scaled_steepest )
            for (j = 0; j < n; j++)
                dx[j] /= anorm[j] * anorm[j];

        // Load Newton step on free variables into dx.
        bcls_scatter( nFree, ix, dxFree, dx );

        // Projected search along dx.
        // Note that ix will be reset.
        if (ls->proj_search == BCLS_PROJ_SEARCH_EXACT)
            err = bcls_proj_search( ls, m, n, x, dx, f, g, aBreak, iBreak,
                                    ix, bl, bu,
                                    ls->wrk_u, ls->wrk_v, ls->wrk_w,
                                    &nSteps );
        else
            err = bcls_proj_backtrack( ls, m, n, x, dx, f, g, bl, bu,
                                       ls->wrk_v, ix, ls->wrk_u,
                                       &nSteps );

        // -------------------------------------------------------------
        // Check for linesearch failures.
        // -------------------------------------------------------------

        // Exit if linesearch found a direction of unbounded descent.
        if (err == BCLS_EXIT_UNBND) {
            ls->exit = err;
            break;
        }

        // Generic linesearch failure.  Further tests needed.
        else if (err == BCLS_EXIT_LFAIL) {

            // If using approx linesearch, revert to exact linesearch.
            if (ls->proj_search == BCLS_PROJ_SEARCH_APPROX) {
                ls->proj_search =  BCLS_PROJ_SEARCH_EXACT;
                PRINT1(" W: Too many approx backtrack iterations."
                       " Switching to exact linesearch.\n");
            }
            // Already using exact linesearch.  Exit.
            else {
                ls->exit = err;
                break;
            }
        }

    } // end of while (1) -- the main iteration loop.

    // Collect various solution statistics.
    ls->soln_obj   = f;
    ls->soln_rNorm = rNorm;
    ls->soln_jInf  = jInf;
    ls->soln_dInf  = dInf;

    // Exit.
    return;
}


/*!

  \brief Compute a Newton step on the free variables.
  
  Compute a Newton step on the variables indexed by ix.  The integer
  array ix (of length nFree) stores the indices of free variables.
  Thus, the kth free variable, with 0 <= k < nFree, is actually
  variable ix[k].

  Conceptually, suppose that x and A are partitioned as
 
      x = [ x_B  x_N ]   and   A = [ B  N ]

  with x_N denoting the free variables.  The objective function is
  given by

      f(x) = 1/2 || Ax - b ||^2 + c'x + 1/2 damp^2 || x ||^2.

  With respect to the free variables, the gradient and Hessian are
  
       g_N(x) = N'(Ax - b) + cF + damp^2 xF = -N'r + cF + damp^2 xF
  and  H_N(x) = N'N + damp^2 I,

  where
  -      r = b - Ax      is the usual residual vector.
  -     cF = c(free)     is the vector of linear terms indexed by ix.
  -     xF = x(free)     is the vector of variables    indexed by ix.
  -      I = eye(nFree).

  The (regularized) Newton step dx_N is then the solution of

                 H_N(x) dx_N = - g_N(x)                       (*)

  or equivalently,

       (N'N + damp^2 I) dx_N = N'r - cF - damp^2 xF,          (**)

  which is really the normal equations for the least-squares problem

         || (      N )       (           r           ) ||
     min || (        ) dx -  (                       ) ||.   (***)
         || ( damp I )       ( - 1/damp cF - damp xF ) ||

  Recall that  I  has dimension (nFree,nFree).

  Inexact Newton
  --------------
  An inexact Newton step is based on approximately satisfying (*).
  The optimality criteria for this step is given by

        || H_N(x) dx_N + g_N(x) || <= rTol || g_N(x) ||,     (****)

  where rTol < 1.  This routine (bcls_newton_step) computes the RHS
  of (****) and passes that to the subroutines as a required
  optimality.

  SUBROUTINES
  -----------
  This routine uses either one of the following
     -  CGLS (via bcls_newton_step_cgls) to solve (**)
     -  LSQR (via bcls_newton_step_lsqr) to solve (***).

  \param[in,out]  ls     The BCLS problem context.
  \param[in]      m      Number of rows in A.
  \param[in]      nFree  Number of columns in N.
  \param[in]      ix     Index of free variables (cols of N).
  \param[in]      damp   Sqrt of regularization parameter on x.
  \param[in]      itnLim Maximum number of minor iterations allowed.
  \param[in]      tol    Required subproblem optimality.
  \param[out]     dxFree The Newton step (solution of the LS subprob).
  \param[in]      x      The current iterate.
  \param[in]      c      The linear term.
  \param[in,out]  r      The current residual.
  \param[out]     itns   No. of iterations taken by the subprob solver.
  \param[out]     opt    The optimality of the subproblem solution.

  \return
  - 0 - Required accuracy was achieved.
  - 1 - The iteration limit was reached.
  - 2 - The matrix is excessively ill-conditioned.
  - 3 - Instability seems likely.

*/
static int
bcls_newton_step( BCLS *ls, int m, int nFree, int ix[], double damp,
		  int itnLim, double tol, double dxFree[], double x[],
		  double c[], double r[], int *itns, double *opt )
{
    int err;
    const int preconditioning = ls->Usolve != NULL;
    double *dx = ls->dx;  // Used as workspace.

    // -----------------------------------------------------------------
    // Initialize the Aprod routine, and the Usolve routine if given.
    // -----------------------------------------------------------------
    bcls_aprod ( ls, BCLS_PROD_INIT, nFree, ix, NULL, NULL );
    if (preconditioning)
        bcls_usolve( ls, BCLS_PRECON_INIT, nFree, ix, NULL, NULL );

    // -----------------------------------------------------------------
    // Print dividing line in the subproblem log.
    // -----------------------------------------------------------------
    if (ls->minor_file) {
        fprintf(ls->minor_file,
                "\n"
                "------------------------- "
                "BCLS Major Iteration: %6d "
                "------------------------- "
                "\n\n", ls->itnMaj);
        fflush(ls->minor_file);
    }
   
    // -----------------------------------------------------------------
    // Call one of the subproblem solvers.
    // -----------------------------------------------------------------
    if (ls->newton_step == BCLS_NEWTON_STEP_LSQR)

        err = bcls_newton_step_lsqr( ls, m, nFree, ix, damp,
                                     itnLim, tol, dxFree, x,
                                     c, r, itns, opt );
    
    else

        err = bcls_newton_step_cgls( ls, m, nFree, ix, damp,
                                     itnLim, tol, dxFree, x,
                                     c, r, itns, opt );

    // -----------------------------------------------------------------
    // If the user is preconditioning, then dxFree currently lives in
    // the preconditioned space.  We need one more call to Usolve to
    // recover dxFree in the original space.  Take this opportunity to
    // call Usolve in its termination mode.
    // -----------------------------------------------------------------
    if (preconditioning) {

        bcls_usolve( ls, BCLS_PRECON_U, nFree, ix, dx, dxFree );
        cblas_dcopy( nFree, dx, 1, dxFree, 1 );  //  dxFree <- dx

        bcls_usolve( ls, BCLS_PRECON_TERM, nFree, ix, NULL, NULL );
    }

    // Call Aprod in its termination mode.
    bcls_aprod ( ls, BCLS_PROD_TERM, nFree, ix, NULL, NULL );

    // Exit.
    return err;
}

/*!

  \brief Find a minimizer along the projected search direction.

  Find the exact minimizer of the quadratic function

    q(t) = 1/2 || A X(t) - b ||^2 + 1/2 damp^2 || X(t) ||^2,

  where X(t) is the projection of (x0 + t dx) onto the box defined by
  bl and bu.
 
  The arguments  y  and  e  are workspace vectors.  They must be at
  least length n and their contents will be overwritten.
 
  \return
  -                0 - No errors encountered.
  - #BCLS_EXIT_UNBND - Encountered direction of infinite descent.

*/
static int
bcls_proj_search( BCLS *ls, int m, int n, double x[], double dx[],
		  double f, double g[], double aBreak[], int iBreak[],
		  int ix[], double bl[], double bu[],
		  double ex[], double Ad[], double Ae[], int *totHits )
{
    int j, k;
    int err     = 0;  // Error flag.
    int nBrkPts = 0;  // No. of breakpoints encountered.
    int nFree   = 0;  // No. of variables that are free to move.
    int xlower;
    int xupper;
    int xfree;
    int stepstat = 0; // Flag if step is
    const int STEP_UND = -1; //   undefined
    const int STEP_LHS =  0; //   on LHS
    const int STEP_MID =  1; //   to interior of interval
    const int STEP_RHS =  2; //   to RHS
    const int STEP_INF =  3; //   unbounded descent.
    int gdpos;        // Flag if gd is "positive".
    int gdzero;       // Flag if gd is "zero".
    int dHdNonNeg;    // Flag if dHd is "nonnegative".
    int dHdNonPos;    // Flag if dHd is "nonpositive".
    int jHits;        // No. of vars that hit ther bound at break point j.

    double lhs;       // The start (LHS) of the current interval.
    double rhs;       // The end (RHS) of the current interval.
    double step;      // Fraction of the current search interval.
    double maxStep;   // The width of the current search interval.
    double q;         // The piecewise-quadratic objective function.
    double gd;        // Projected gradient of q along current arc.
    double ge;        // Projected gradient of q along last arc.
    double dHd;       // Projected Hessian  of q along current arc.
    double dHe;       // Inner product of Hessian along d and update e.

    const int damped = ls->damp > 0.0;
    const double damp2 = ls->damp * ls->damp;
    const double eps = ls->eps;
    const double epsx = ls->epsx;
    const double epsfixed = ls->epsfixed;
    const double BigNum = ls->BigNum;

    PRINT6("\n%3s %10s %10s %10s %10s %10s\n",
	   "Var","Lower","x","Upper","dx","Break");

    // -----------------------------------------------------------------
    // Find out which variables are free to move.
    // -----------------------------------------------------------------
    for (j = 0; j < n; j++) {

	if ( bu[j] - bl[j] < epsfixed ) {
	    // This variable is fixed.  Zero out its search direction.
	    dx[j] = 0.0;
	    continue;
	}
	
	xupper = bu[j] - x[j] <= epsx; // Variable on upper bound.
	xlower = x[j] - bl[j] <= epsx; // Variable on lower bound.
	xfree  = !xlower && !xupper;   // Variable is floating btwn bnds.
	
	assert( xupper || xlower || xfree );

	if ( xfree ) {
	    if ( fabs(dx[j]) > eps ) {
		iBreak[nFree] = j;
		nFree++;
	    }
	    else
		dx[j] = 0.0;
	}

	else if ( xupper ) {
	    if ( dx[j] < - eps ) {
		iBreak[nFree] = j;
		nFree++;
	    }
	    else
		dx[j] = 0.0;
	}

	else { // if ( xlower )
	    if ( dx[j] > eps ) {
		iBreak[nFree] = j;
		nFree++;
	    }
	    else
		dx[j] = 0.0;
	}
    }
    
    // -----------------------------------------------------------------
    // Compute the breakpoints of the free variables.
    //
    // Each breakpoint measures the maximum allowable step that will
    // take the corresponding free variable to its boundary.
    // -----------------------------------------------------------------
    for (j = 0; j < nFree; j++) {
	k = iBreak[j];
	if ( dx[k] > 0.0 )
	    aBreak[j] = ( bu[k] - x[k] ) / dx[k];
	else
	    aBreak[j] = ( bl[k] - x[k] ) / dx[k];
    }

    // -----------------------------------------------------------------
    // Diagnostic printing.
    // -----------------------------------------------------------------
    if (ls->print_level >= 6) {
	k = 0;
	for (j = 0; j < n; j++) {
	    PRINT6("%3d %10.3e %10.3e %10.3e %10.3e ",
		   j, bl[j], x[j], bu[j], dx[j]);
	    if ( j == iBreak[k] ) {
		PRINT6("%10.3e\n",aBreak[k]);
		k++;
	    }
	    else
		PRINT6("%10s\n","--");
	}
	PRINT6("\n");
    }

    // Exit if no variables are free to move.
    if ( !nFree ) return 0;

    // No. of breakpoints remaining is initially the number of free vars.
    nBrkPts = nFree;

    // Build the breakpoints into a heap.
    bcls_heap_build( nBrkPts, aBreak, iBreak );

    // -----------------------------------------------------------------
    // Projected gradient of q along dx.
    // -----------------------------------------------------------------
    gd = 0.0;
    for (j = 0; j < nFree; j++) {
	k = iBreak[j];
	gd += g[k] * dx[k];
    }

    // Quick exit if the directional derivative is positive.
    if ( gd > epsx ) return 0;

    // -----------------------------------------------------------------
    // Projected Hessian of q along dx:
    // dHd = dx' A' A dx = (A dx)' (A dx) = Ad' Ad.
    // -----------------------------------------------------------------
    bcls_aprod( ls, BCLS_PROD_A, nFree, iBreak, dx, Ad );
    dHd = cblas_ddot( m, Ad, 1, Ad, 1 );
    if (damped)
	dHd += damp2 * cblas_ddot( n, dx, 1, dx, 1 );

    // Initialization.
    q        = 0.0; // Value of piecewise quadratic.
    rhs      = 0.0; // This initializes the very first breakpoint.
    *totHits = 0;   // Initially no hits.

    // -----------------------------------------------------------------
    // Explore each piecewise quadratic intervals until non remain.
    // -----------------------------------------------------------------
    while (nBrkPts > 0) {

	// Set the status of the current step to "undefined".
	stepstat = STEP_UND;

	// Grab the boundaries of the current search interval.
	lhs     = rhs;
	rhs     = aBreak[0];  // The next breakpoint (top of the heap).
	maxStep = rhs - lhs;  // The maximum step.
	k       = iBreak[0];  // Var corresponding to the current break.

	// The smallest breakpoint has been sampled.  Move it to the
	// top of the heap and restore the heap property.
	nBrkPts = bcls_heap_del_min( nBrkPts, aBreak, iBreak );

	// -------------------------------------------------------------
	// Test the various exit or continuation conditions.
	// -------------------------------------------------------------
	gdpos     = gd       >    epsx;
	gdzero    = fabs(gd) <=   epsx;
	dHdNonNeg = dHd      >= - epsx;
	dHdNonPos = dHd      <    epsx;

	if ( gdpos || ( gdzero && dHdNonNeg) ) {
	    // ---------------------------------------------------------
	    // Minimimum is on LHS. (Captures the  dx = 0  case.)  Done.
	    // ---------------------------------------------------------
	    step = 0.0;
	    stepstat = STEP_LHS;
	    break;
	}	

	// If we got this far, gd is zero or negative.

	else if ( dHdNonPos ) {
	    // ---------------------------------------------------------
	    // Found dir of neg curvature.  Minimizer is either on RHS,
	    // or at -infinity.  Check.
	    // ---------------------------------------------------------
	    if ( maxStep >= BigNum ) {
		stepstat = STEP_INF;
		break;
	    }
	    else {
		step = maxStep;
		stepstat = STEP_RHS;
	    }
	}
	else {
	    // ---------------------------------------------------------
	    // Minimum is on RHS or inside the interval.  Check.
	    // ---------------------------------------------------------
	    step = - gd / dHd;
	    if ( step >= BigNum  &&  step < maxStep ) {
		// Step is unbounded.  Exit with unbounded flag.
		stepstat = STEP_INF;
		break;
	    }
	    else if ( step <  maxStep ) {
		// Minimum found in the interval's interior.  Done.
		stepstat = STEP_MID;
		break;
	    }
	    else if ( step >= maxStep ) {
		// Minimum is on RHS.  Truncate the step and continue.
		step = maxStep;
		stepstat = STEP_RHS;
	    }
	}
	// -------------------------------------------------------------
	// Take the step.
	// -------------------------------------------------------------
	cblas_daxpy( n, step, dx, 1, x, 1 );
	    
	// -------------------------------------------------------------
	// Other variables may share (essentially) the same breakpoint.
	//
	// Loop over these (and the current) variables, pushing them
	// onto their bounds, zeroing out the corresponding search
	// direction.
	// -------------------------------------------------------------
	jHits = 0; // No. of vars that hit their bound at this breakpoint.
	while(1) {
	    
	    // There may have been roundoff error.  Push the var that
	    // hit its breakpoint exactly onto its bound.
	    if (dx[k] < 0.0)
		x[k] = bl[k];
	    else
		x[k] = bu[k];

	    // Record the index and step of the var that just hit its bound.
	    ex[k] = dx[k];
	    ix[jHits] = k;
	    jHits++;

	    // That var no longer contributes.  Zero out its step.
	    dx[k] = 0.0;
	    
	    // Print something out.
	    PRINT6("%3d %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
		   k,bl[k],x[k],bu[k],dx[k],gd,dHd);

	    // No more breakpoints left to explore.
	    if (nBrkPts == 0)
		break;

	    // Next breakpoint is significantly different.
	    else if (aBreak[0] > rhs + 10*eps)
		break;

	    // Next breakpoint is essentially the same.
	    else {
		k = iBreak[0];
		nBrkPts = bcls_heap_del_min( nBrkPts, aBreak, iBreak );
	    }
	}
	
	// Update the aggregate no. of breakpoints encountered so far.
	*totHits += jHits; 

	// -------------------------------------------------------------
	// Update the function value.  (Used for printing only.)
	// -------------------------------------------------------------
	q  = q + step * (gd + 0.5 * step * dHd);

	// -------------------------------------------------------------
	// Update the directional derivative and Hessian for the next
	// search direction.
	// -------------------------------------------------------------
	
	// Inner products: g'e.
	ge = 0.0;
	for (j = 0; j < jHits; j++) {
	    k = ix[j];
	    ge += g[k]*ex[k];
	}

	// Inner product: d' H e.
	bcls_aprod( ls, BCLS_PROD_A, jHits, ix, ex, Ae ); //  Ae = A * ex
	dHe = cblas_ddot( m, Ad, 1, Ae, 1 );              // dHe = d' A' A e
	if (damped) {
	    for (j = 0; j < jHits; j++) {
		k = ix[j];
		dHe += damp2 * dx[k] * ex[k];
	    }
	}

	// Update projected gradient.
	gd  = gd - ge + step * ( dHd - dHe );
	
	// Update projected Hessian.
	cblas_daxpy( m, -1.0, Ae, 1, Ad, 1 );     //  Ad <- Ad - Ae
	dHd = cblas_ddot( m, Ad, 1, Ad, 1 );      // dHd <- d' A' A d
	if (damped)
	    dHd = dHd + damp2 * cblas_ddot( n, dx, 1, dx, 1 );

    } // end of while (nBrkPs > 0)

    if ( stepstat == STEP_MID )
	// Step into the middle of an interval and exit.
	cblas_daxpy( n, step, dx, 1, x, 1 );

    else if ( stepstat == STEP_INF )
	err = BCLS_EXIT_UNBND;

    // Return the status of the last step.
    return err;
}

/*!

  \brief Simple backtracking to satisfy an Armijo condition.

  Do a simple backtracking search on the quadratic function

    q(x + s) = f(x) + g(x)'s + 0.5 s' H s,

  where H = A'A + damp^2 I and s = P[dx] is the projection of dx onto
  the box defined by bl and bu.  The function value f = f(x) and the
  gradient g = g(x) are defined at entry and left unchanged.  The step
  dx is reduced by a contant amount at each iteration until its
  projection s satisfies an Armijo condition, ie, sufficient descent is
  achieved when

    g's + 0.5 s' H s <= mu g's,

  and mu is a fixed constant.

  \return
  -                0 - No errors encountered.
  - #BCLS_EXIT_UNBND - Encountered direction of infinite descent.
  - #BCLS_EXIT_LFAIL - Linesearch failure.

*/
static int 
bcls_proj_backtrack( BCLS *ls, int m, int n, double x[], double dx[],
		     double f, double g[], double bl[], double bu[],
		     double s[], int ix[], double As[], int *nSteps )
{
    int j;
    int armijo = 0;   // Flag: Armijo condition satisfied.
    int unbounded;    // Flag: Step may be unbounded.
    int descent;      // Flag: Descent direction.
    int negCurv;      // Flag: Nonpositive curvature.
    int xupper;       // Flag: Variable on upper bound.
    int xlower;       // Flag: Variable on lower bound.
    int xdecr;        // Flag: Variable is decreasing.
    int xincr;        // Flag: Variable is increasing.
    int nBreakPts;    // No. of break points.
    int jBreakMin;    // Index of variable with the nearest breakpoint.
    double Break;     // Step length to the break point for variable j.
    double BreakMin;  // Step length to the nearest breakpoint.
    double q;         // The quadratic objective function.
    double gs;        // Projected gradient along current search direction.
    double sHs;       // Projected Hessian  along current search direction.
    double step = 1;  // Current step length.

    const int    backlimit = ls->backtrack_limit;
    const int    damped    = ls->damp > 0.0;
    const double damp2     = ls->damp * ls->damp;
    const double mu        = ls->mu;
    const double backtrack = ls->backtrack;
    const double BigNum    = ls->BigNum;
    const double eps       = ls->eps;
    const double epsx      = ls->epsx;

    *nSteps = 0;      // No. of backtracking steps.

    // -----------------------------------------------------------------
    // Need to include contribution for *all* columns of A.
    // -----------------------------------------------------------------
    for (j = 0; j < n; j++) ix[j] = j;

    // -----------------------------------------------------------------
    // Find first breakpoint, and store as BreakMin.
    // -----------------------------------------------------------------
    nBreakPts = 0;
    BreakMin  = BigNum;
    for (j = 0; j < n; j++) {
        xupper = bu[j] - x[j] <= epsx; // Variable on upper bound.
        xlower = x[j] - bl[j] <= epsx; // Variable on lower bound.
        xdecr  = dx[j] < -eps;         // Variable is decreasing.
        xincr  = dx[j] >  eps;         // Variable is increasing.

        if ( !xlower  &&  xdecr ) {
            nBreakPts++;
            Break = ( bl[j] - x[j] ) / dx[j];
            if ( Break < BreakMin ) {
                BreakMin  = Break;
                jBreakMin = j;
            }
        }
        else if ( !xupper  &&  xincr ) {
            nBreakPts++;
            Break = ( bu[j] - x[j] ) / dx[j];
            if ( Break < BreakMin ) {
                BreakMin  = Break;
                jBreakMin = j;
            }
        }
    }
    
    if ( nBreakPts == 0 )
        BreakMin = 0.0;

    PRINT3("I: Proj Backtrack: Steplength to the nearest bound: %9.2e\n",
           BreakMin);
    
    // -----------------------------------------------------------------
    // Backtracking loop.
    // -----------------------------------------------------------------
    while ( !armijo  &&  step > BreakMin ) {

	if (*nSteps >= backlimit)
	    return BCLS_EXIT_LFAIL;
	
        // -------------------------------------------------------------
        // Store the projected step into s:  s = P[x + step*dx] - x.
        // -------------------------------------------------------------
        bcls_project_step( n, s, step, x, dx, bl, bu );

	// -------------------------------------------------------------
	// Compute q(s) = g's + 0.5 s' H s.
	// -------------------------------------------------------------
	bcls_aprod( ls, BCLS_PROD_A, n, ix, s, As ); //  As = A * s
	sHs = cblas_ddot( m, As, 1, As, 1 );         // sHs = (A s)'(A s)
	if (damped)
	    sHs = sHs + damp2 * cblas_ddot( n, s, 1, s, 1 );
	gs  = cblas_ddot( n, g, 1, s, 1 );    //  gs = g's
	q   = gs + 0.5 * sHs;                 //   q = g's + 0.5 s' H s
	
	// -------------------------------------------------------------
	// If sHs <= 0 and gs < 0, there may be a direction of
	// unbounded descent.  Check.
	// -------------------------------------------------------------
	descent  =  gs  <= -epsx;
	negCurv  =  sHs <=  epsx;

	if ( descent  &&  negCurv ) {

	    unbounded = 1;
	    for (j = 0; j < n; j++) {
		if (s[j] > eps)
		    unbounded = bu[j] >=   BigNum;
		
		else if (s[j] < eps)
		    unbounded = bl[j] <= - BigNum;
		
		if (!unbounded)
                    break;
	    }
	    
	    if (unbounded)
                return BCLS_EXIT_UNBND;
	}

	// -------------------------------------------------------------
	// Check if the sufficient decrease criteria was satisfied.
	// Otherwise, reduce the step and continue.
	// -------------------------------------------------------------
        PRINT3("I: Proj Backtrack: step = %9.2e  q = %10.2e  gs = %10.2e\n",
               step, q, gs);

	armijo = ( q <= mu * gs );

	if (!armijo) {
	    step = backtrack * step;
	    *nSteps += 1;
	}
    }

    // -----------------------------------------------------------------
    // Force a variable onto its nearest bound: If the armijo
    // condition is satisfied by a truncated step that falls short of
    // the nearest breakpoint, then we need to be a bit more
    // aggressive and continue the step until the nearest boundary.
    //
    // 17 Dec 06: Removed this step.  Why was I doing it?!
    // -----------------------------------------------------------------
#ifdef UNDEF
    if ( step < 1.0  &&  step < BreakMin ) {
        step = BreakMin;
        PRINT3("I: Proj Backtrack: Force index %4i onto bound\n",
               jBreakMin );
    }
#endif

    // Compute the final projected step.
    bcls_project_step( n, s, step, x, dx, bl, bu );
    
    // Take the step: x = x + s.
    cblas_daxpy( n, 1.0, s, 1, x, 1 );
    
    // Push x onto bounds that are very close.
    for (j = 0; j < n; j++) {
        if      (x[j] < bl[j] + 10*eps)
                 x[j] = bl[j];
        else if (x[j] > bu[j] - 10*eps)
                 x[j] = bu[j];
    }

    return 0;
}
