/* bccgls.c
   $Revision: 281 $ $Date: 2006-12-13 09:03:53 -0800 (Wed, 13 Dec 2006) $

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
  This file implements conjugate gradients for the normal equations
  (CGLS).
*/

#include <cblas.h>
#include <string.h>
#include <assert.h>
#include <float.h>

#include "bcls.h"
#include "bclib.h"
#include "bccgls.h"

/*!

  \brief Mat-vec routine called by CGLS.
  
  CGLS call this routine, which in turn calls bcls_aprod:
  - If      mode == BCLS_PROD_A,   y      = A dxFree.
  - else if mode == BCLS_PROD_At,  dxFree = A'y.

  This routine is declared "static" so that it won't be confused with
  the user's own Aprod routine.
  
  \param[in,out] ls      BCLS solver context.
  \param[in]     mode    Specifies if the product with A or A' is required.
  \param[in]     m       Number of rows in A.
  \param[in]     nFree   Length of dxFree (and ix which is recovered via ls).
  \param[in,out] dxFree  RHS or LHS (depending on the mode).
  \param[in,out] y       RHS or LHS (depending on the mode).

*/
static void
aprod_free( BCLS *ls, const int mode, const int m, const int nFree,
	    double dxFree[], double y[] )
{
    int    *ix = ls->ix;
    double *dx = ls->dx;  // Used as workspace.

    assert( mode == BCLS_PROD_A || mode == BCLS_PROD_At );

    if (mode == BCLS_PROD_A) {

        bcls_scatter( nFree, ix, dxFree, dx );          // dx(ix) = dxFree.
	bcls_aprod( ls, BCLS_PROD_A, nFree, ix, dx, y );// y = A dx.

    } else {
	
	bcls_aprod( ls, BCLS_PROD_At, nFree, ix, dx, y );// dx = A' y.
        bcls_gather( nFree, ix, dx, dxFree );            // dxFree = dx(ix);
    }

    return;
}

/*!

  \brief CG for least squares.

  This routine implements CG for solving the symmetric linear system

     ( N'N + damp I ) x = N'r + c

  where N is an m-by-n submatrix of the general matrix A.  This
  implementation is adapted from Hestenes and Stiefel's CG for
  least-squares method to accomodate the damping term and the
  additional vector on the RHS.  (See, eg, Bjork, Numerical Methods
  for Least-Squares Problems, 1996, p. 289.)

  \param[in,out] ls     BCLS solver context.
  \param[in]     m      No. of rows in A.
  \param[in]     n      No. of cols in A.
  \param[in]     kmax   Maximum no. of CG iterations.
  \param[in]     tol    Required optimality tolerance:
                        - norm( A'A x - A'b ) / norm(A'b)
  \param[in]     damp   Regularization term on x.
  \param[in]     c      Additional vector (set NULL if not present).
  \param[out]    x      Solution vector.
  \param[in,out] r      On entry, r must be defined.\n
                        On exit,  r is the residual.
  \param[in,out] s      Work vector.
  \param[in,out] p      Work vector.
  \param[out]    itns   No. of CG iterations.
  \param[out]    opt    Relative residual (see tol, above).
  \param[in]     nout   File for log output (set NULL if not needed).

  \returns
  - 0:  Required accuracy was achieved.
  - 1:  The iteration limit (kmax) was reached.
  - 2:  N'N appears to be singular.
  
  \todo Update normx rather than recompute at each iteration.

*/
static int
bcls_cgls( BCLS *ls, int m, int n, int kmax, double tol, double damp,
	   double c[], double x[], double r[], double s[], double p[],
	   int *itns, double *opt, FILE *nout )
{    
    int
        k = 0, info = 0, converged = 0, unstable = 0;
    double
        xmax = 0.0, resNE, Arnorm, Arnorm0,
        normp, norms, normx = 0.0, alpha, beta,
        gamma, gamma1, delta;
    const int
        damped = damp > 0.0;
    const double
        eps = DBL_EPSILON;

    // Initialize.
    bcls_dload( n, 0.0, x, 1 );                  // x  = 0
    aprod_free( ls, BCLS_PROD_At, m, n, s, r );  // s  = A'r
    if (c) 
        cblas_daxpy( n, 1.0, c, 1, s, 1 );       // s += c
    cblas_dcopy( n, s, 1, p, 1 );                // p  = s
    Arnorm0    = cblas_dnrm2( n, s, 1 );
    gamma      = Arnorm0 * Arnorm0;

    if (nout) {
        fprintf(nout,"\n CGLS:\n");
        fprintf(nout,"  Rows............. %8d\t", m);
        fprintf(nout,"  Columns.......... %8d\n", n);
        fprintf(nout,"  Optimality tol... %8.1e\t", tol);
        fprintf(nout,"  Iteration limit.. %8d\n", kmax);
        fprintf(nout,"  Damp............. %8.1e\t", damp);
        fprintf(nout,"  Linear term...... %8s\n",
                (c) ? "yes" : "no" );
        fprintf(nout, "\n\n%5s %17s %17s %9s %10s\n",
                "k", "x(1)","x(n)","normx","resNE");
        fprintf(nout, "%5d %17.10e %17.10e %9.2e %10.2e\n",
                k, x[0], x[n-1], normx, 1.0 );
    }
    // -----------------------------------------------------------------
    // Main loop.
    // -----------------------------------------------------------------
    while ( k < kmax  &&  !converged  &&  !unstable ) {

        k++;

        aprod_free( ls, BCLS_PROD_A, m, n, p, s );   // s = A p

        norms = cblas_dnrm2( m, s, 1 );
        delta = norms * norms;
        if (damped) {
            normp  = cblas_dnrm2( n, p, 1 );
            delta += damp * normp * normp;
        }
        if (delta <= eps)
            delta =  eps;

        alpha = gamma / delta;

        cblas_daxpy( n,  alpha, p, 1, x, 1 );        // x  = x + alpha*p
        cblas_daxpy( m, -alpha, s, 1, r, 1 );        // r  = r - alpha*s
        aprod_free( ls, BCLS_PROD_At, m, n, s, r );  // s  = A'r
        if (damped)
            cblas_daxpy( n, -damp, x, 1, s, 1 );     // s += - damp*x
        if (c)
            cblas_daxpy( n,   1.0, c, 1, s, 1 );     // s += c

        Arnorm = cblas_dnrm2( n, s, 1 );
        gamma1 = gamma;
        gamma  = Arnorm * Arnorm;
        beta   = gamma / gamma1;
        
        cblas_dscal( n, beta, p, 1 );                // p  = beta*p
        cblas_daxpy( n, 1.0, s, 1, p, 1);            // p += s 

        // Check convergence.
        normx     = cblas_dnrm2( n, x, 1 );
        xmax      = fmax( xmax, normx );
        converged = Arnorm <= Arnorm0 * tol;
        unstable  = normx * tol >= 1.0;

        // Output.
        resNE = Arnorm / Arnorm0;
        if (nout)
            fprintf(nout, "%5d %17.10e %17.10e %9.2e %10.2e\n",
                    k, x[0], x[n-1], normx, resNE );
        
    } // while

    double shrink = normx / xmax;
    if (converged)
        info = 0;
    else if (k >= kmax)
        info = 1;
    else if (shrink <= sqrt(tol))
        info = 2;

    // Output.
    if (nout) fprintf(nout,"\n");
    *itns = k;
    *opt  = resNE;
    return info;
}

/*!

  \brief Compute a Newton step using CGLS.

  See the description of #bcls_newton_step.

  \param[in,out] ls      BCLS solver context
  \param[in]     m       No. of rows in A.
  \param[in]     nFree   No. of free columns, i.e., no. of cols in N.
  \param[in]     ix      Index of vars for which a Newton step is computed.
  \param[in]     damp    Regularization term on dxFree.
  \param[in]     itnLim  Maximum no. of CG iterations.
  \param[in]     tol     Required optimality tolerance:
                         - norm( N'N x - N'b ) / norm(N'b).
  \param[in]     dxFree  Newton direction in the free vars (len = nFree).
  \param[in]     x       Current iterate.
  \param[in]     c       Linear term.  (Set to NULL if not present.)
  \param[in,out] r       On entry, r is current residual: r = b - Ax. \n
                         Will be used as workspace.
  \param[in]     itns    No. of CG iterations.
  \param[in]     opt     Relative residual (see tol, above).

  \return The termination flag returned by #bcls_cgls.

*/
int
bcls_newton_step_cgls( BCLS *ls, int m, int nFree, int ix[], double damp,
                       int itnLim, double tol, double dxFree[], double x[],
                       double c[], double r[], int *itns, double *opt )
{
    int
        j;
    double
        *u;
    const double
        damp2 = damp*damp;

    if (c) {
        u = ls->wrk_u;
        bcls_gather( nFree, ix, c, u );       // u  =   c(ix).
        cblas_dscal( nFree, -1.0, u, 1 );     // u  = - c(ix).
        if ( damp > 0.0 )                     // u += - damp2 x(ix).
            for (j = 0; j < nFree; j++)
                u[j] -= damp2 * x[ix[j]];
    }
    else
        u = NULL;

    return bcls_cgls( ls, m, nFree, itnLim, tol, damp2,
                      u, dxFree, r, ls->wrk_v, ls->wrk_w,
                      itns, opt, ls->minor_file );
}
