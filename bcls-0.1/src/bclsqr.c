/* bclsqr.c
   $Revision: 273 $ $Date: 2006-09-04 15:59:04 -0700 (Mon, 04 Sep 2006) $

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
   Interface to LSQR routine.  Used to compute a Newton step.
*/

#include <cblas.h>
#include <string.h>
#include <stdio.h>

#include "bcls.h"
#include "bclib.h"
#include "bclsqr.h"
#include "lsqr.h"

/*!

  \brief Mat-vec routine called by LSQR.

  LSQR calls this routine, which in turn calls bcls_aprod: This
  routine is declared "static" so that it won't be confused with the
  user's own Aprod routine.
  
  - If mode = 1,
    - y(1:m)  <-  y(1:m)  +  A(:,ix) * dxFree.
    - y(m+1:) <-  y(m+1:) +     damp * dxFree

  - If mode = 2,
    - dxFree <- dxFree + A'* y(1:m) + damp y(m+1:).

  Note that mSubProb is the number of rows in the subproblem.  This
  may or may not be equal to m, which is the number of rows in the
  original problem.  If the subproblem is damped (because the user has
  either provided a linear term or an explicit damping parameter),
  then mSubProb = m + nFree.  Otherwise, mSubProb = m.

  \param[in]      mode      Determines which producte with A is required.

  \param[in]      mSubProb  Number of rows in the matrix seen by LSQR.  Also:
                            - length of y.
  \param[in]      nFree     Number of columns in A(:,ix).  Also:
                            - length of ix
                            - length of dxFree.
  \param[in,out]  dxFree    Primal variables
  \param[in,out]  y         Dual variables.
  \param[in,out]  UsrWrk    Transit pointer to the BCLS problem context.

*/
static void
aprod_free_lsqr( const int mode, const int mSubProb, const int nFree,
                 double dxFree[], double y[], void *UsrWrk)
{
    int j;
    BCLS *ls = (BCLS *)UsrWrk;  // Reclaim access to the BCLS workspace.
    const int    m      = ls->m;
    const int    preconditioned = ls->Usolve != NULL;
    const int    damped = ls->damp_actual > 0.0;
    const double damp   = ls->damp_actual;
    int          *ix    = ls->ix;
    double       *dx    = ls->dx;    // Used as workspace.
    double       *dy    = ls->wrk_u; // ...

    if (mode == 1) {
        
        // Solve U dy = dxFree.  U is nFree-by-nFree, and so the
        // relevant part of dy has length nFree.
        if (preconditioned)
            bcls_usolve( ls, BCLS_PRECON_U, nFree, ix, dy, dxFree );
        else
            cblas_dcopy( nFree, dxFree, 1, dy, 1 );

        // y2 <- y2 + damp * dy(1:nFree).
        if (damped)
            cblas_daxpy( nFree, damp, dy, 1, &y[m], 1 );
        
        // Scatter dy into dx: dx(ix) <- dy(1:nFree).
        for (j = 0; j < nFree; j++)  dx[ ix[j] ] = dy[j];
        
        // dy <- A(:,ix) * dx(ix).
        bcls_aprod( ls, BCLS_PROD_A, nFree, ix, dx, dy );
        
        // y1 <- y1 + dy.
        cblas_daxpy( m, 1.0, dy, 1, y, 1 );
        
    }
    else { // mode == 2
        
        // dx <- A' * y1.
        bcls_aprod( ls, BCLS_PROD_At, nFree, ix, dx, y );
        
        // Gather dx(ix) into dy: dy(1:nFree) <- dx(ix).
        for (j = 0; j < nFree; j++) dy[j] = dx[ ix[j] ];
        
        // dy <- dy + damp * y2.
        if (damped)
            cblas_daxpy( nFree, damp, &y[m], 1, dy, 1 );
        
        // Solve U' dx = dy.
        if (preconditioned)
            bcls_usolve( ls, BCLS_PRECON_Ut, nFree, ix, dy, dx );
        else
            cblas_dcopy( nFree, dy, 1, dx, 1 );

        // dxFree <- dxFree + dx.
        cblas_daxpy( nFree, 1.0, dx, 1, dxFree, 1 );
        
    }
    return;
}

/*!

  \brief Compute a Newton step using LSQR.

  \see  bcls_newton_step_cgls.

  \param[in,out] ls      BCLS problem context.
  \param[in]     m       Number of rows in A.
  \param[in]     nFree   Number of columns in A(:,ix).  Also:
                         - length of ix and dxFree
  \param[in]     ix      Index of free variables.
  \param[in]     damp    Regularization parameter.
  \param[in]     itnLim  Iteration limit on current LSQR call.
  \param[in]     tol     LSQR's atol and btol.
  \param[in,out] dxFree  Search direction on free variables.
  \param[in,out] x       Current point.
  \param[in]     c       Linear term.
  \param[in,out] r       Residual. Used as RHS for LSQR.
  \param[out]    itns    Number of LSQR iterations on current subproblem.
  \param[out]    opt     Optimality achieved by LSQR on current subproblem.

  \return
  - 0: Required accurace was achieved.
  - 1: The iteration limit (itnLim) was reached.
  - 2: A(:,ix) is excessively ill-conditioned.

*/
int
bcls_newton_step_lsqr( BCLS *ls, int m, int nFree, int ix[], double damp,
		       int itnLim, double tol, double dxFree[], double x[],
		       double c[], double r[], int *itns, double *opt )
{
    int j, k;                             // Misc. counters.
    int mpn;                              // No. of rows in [ N; damp I ].
    int          unscale_dxFree   = 0;
    const int    rescaling_method = 1;
    const int    linear   = c != NULL;
    const int    damped   = damp  > 0.0;
    const double damp_min = ls->damp_min;
    const double damp2    = damp * damp;
    const double zero     = 0.0;
    double damp_actual;
    double beta1, beta2;

    // LSQR outputs
    int    istop;      // Termination flag.
    double anorm;      // Estimate of Frobenious norm of Abar.
    double acond;      // Estimate of condition no. of Abar.
    double rnorm;      // Estimate of the final value of norm(rbar).
    double xnorm;      // Estimate of the norm of the final solution dx.

    // Set r(m+1:) <- - 1/beta c + damp^2/beta x, where
    //     beta = max(min_damp, damp).
    // Note that r is declared length (m+n), so there is always enough
    // space.  The chosen damping parameter is stored in
    // ls->damp_actual so that it can be used in bcls_aprod_free.
    if (!damped  &&  !linear) {
	mpn = m;
	ls->damp_actual = 0.0;
    }
    else {

	mpn = m + nFree;

	//--------------------------------------------------------------
	// Rescale the RHS.  Will need to unscale LSQR's solution later.
	//--------------------------------------------------------------
	if (rescaling_method) {

	    if (linear) {

		unscale_dxFree  = 1;
		damp_actual     = fmax( damp, damp_min );
		ls->damp_actual = damp_actual;
		
		// r(1:m) = damp_actual * r(1:m)
		cblas_dscal( m, damp_actual, r, 1 );

		// r(m+1:) = - c
		for (j = 0; j < nFree; j++ )
		    r[m+j] = - c[ ix[j] ];
		
		// r(m+1:) = - c - damp^2 * x
		if (damped)
		    for (j = 0; j < nFree; j++)
			r[m+j] -= damp2 * x[ ix[j] ];
	    }
	    else {
		
	        ls->damp_actual = damp;
		for (j = 0; j < nFree; j++)
		    r[m+j] = - damp * x[ ix[j] ];
	    }
	}
	//--------------------------------------------------------------
	// No rescaling.
	//--------------------------------------------------------------
	else {
	    if (damped  &&  linear) {
		
		damp_actual     =  fmax( damp, damp_min );
		beta1           =  -1.0 / damp_actual;
		beta2           =  -damp * damp / damp_actual;
		ls->damp_actual =  damp_actual;
		
		for (j = 0; j < nFree; j++) {
		    k = ix[j];
		    r[m+j] = beta1 * c[k] + beta2 * x[k];
		}
	    }
	    else if ( damped  &&  !linear ) {
		
		beta2 = - damp;
		ls->damp_actual = damp;
		for (j = 0; j < nFree; j++)
		    r[m+j] = beta2 * x[ ix[j] ];
	    }
	    else if (!damped  &&   linear ) {
		
		beta1 = - 1.0 / damp_min;
		ls->damp_actual = damp_min;
		for (j = 0; j < nFree; j++)
		    r[m+j] = beta1 * c[ ix[j] ];
	    }
	}
    }
    
    // -----------------------------------------------------------------
    // Solve the subproblem with LSQR.
    // -----------------------------------------------------------------
    bcls_timer( &(ls->stopwatch[BCLS_TIMER_LSQR]), BCLS_TIMER_START );

    lsqr( mpn, nFree, aprod_free_lsqr, zero, (void *)ls,
          r, ls->wrk_v, ls->wrk_w, dxFree, NULL,
          tol, tol, ls->conlim, itnLim, ls->minor_file,
          &istop, itns, &anorm, &acond,
          &rnorm, opt, &xnorm );

    bcls_timer( &(ls->stopwatch[BCLS_TIMER_LSQR]), BCLS_TIMER_STOP );

    // -----------------------------------------------------------------
    // Cleanup and exit.  First unscale dxFree if needed.
    // -----------------------------------------------------------------
    if (rescaling_method  &&  unscale_dxFree)
	cblas_dscal( nFree, 1.0/damp_actual, dxFree, 1 );

    if (istop <= 3)
	return 0;
    else if (istop == 5)
	return 1;
    else
	return 2;
}
