/* bclib.c
   $Revision: 285 $ $Date: 2006-12-20 09:33:50 -0800 (Wed, 20 Dec 2006) $

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
   Helper routines for BCLS.
*/

#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <setjmp.h>

#include "bclib.h"

/*!

  \brief Compute a projected step.
  
  Given a current point  x  and a step  dx, compute the projected step
  s = P[x + step*dx] - x.
  
  \param[in]   n      Length of x, s, dx, bl, bu.
  \param[out]  s      The projected step.
  \param[in]   step   Step length.
  \param[in]   x      The current point.
  \param[in]   dx     The step that needs to be projected.
  \param[in]   bl     Lower bounds on x.
  \param[in]   bu     Upper bounds on x.
  
*/
void
bcls_project_step( int n, double s[], double step, double x[],
                   double dx[], double bl[], double bu[] )
{
    int j;
    double xpj;

    for (j = 0; j < n; j++) {

        xpj = x[j] + step*dx[j];

        if ( xpj < bl[j] )
            s[j] = bl[j] - x[j];
        
        else if ( xpj > bu[j] )
            s[j] = bu[j] - x[j];
        
        else
            s[j] = step * dx[j];
    }
    return;
}

/*!

  \brief Examine the bounds bl and bu.

  Print some statistics about the problem bounds, and make sure that
  the problem is feasible, i.e., that bl <= bu.

  \param[in] ls  BCLS solver context.
  \param[in] n   Length of bl and bu.
  \param[in] bl  Lower bounds on x.
  \param[in] bu  Upper bounds on x.

  \return
  -                0: Bounds are   feasible.
  - #BCLS_EXIT_INFEA: Bounds are infeasible.

*/
int
bcls_examine_bnds( BCLS *ls, const int n, double bl[], double bu[] )
{
    int j, jopen;
    double blj, buj;
    const double bigL = - ls->BigNum;
    const double bigU =   ls->BigNum;
    const double epsfixed = ls->epsfixed;

    // Variable-type counts.  No. of variables that
    int neql = 0; // are "open" (have an interior)
    int nfix = 0; // are fixed
    int nneg = 0; // are unbounded below
    int npos = 0; // are unbounded above
    int nlow = 0; // have a finite lower bound
    int nupp = 0; // have a finite upper bound
    int nnrm = 0; // have finite lower and upper bounds
    int nfre = 0; // are free.

    for (j = 0; j < n; j++) {
	jopen = 0;
	blj   = bl[j];
	buj   = bu[j];
	
	if (buj < blj)                  // Infeasible interval.  Exit.
	    goto infeas;
	
	else if (buj - blj <= epsfixed) // Fixed variable.
	    nfix++;

	else {                          // Open interval.
	    neql++;
	    jopen = 1;
	}

	if     (blj <= bigL  &&  buj == 0            ) nneg++;
	if     (blj == 0     &&  buj >= bigU         ) npos++;
	if     (blj >  bigL                  && jopen) nlow++;
	if     (                 buj <  bigU && jopen) nupp++;
	if     (blj >  bigL  &&  buj <  bigU         ) nnrm++;
	if     (blj <= bigL  &&  buj >= bigU         ) nfre++;
    }
    
    ls->unconstrained = nfre == n;

    /*    PRINT1("\n Bounds:");
    if (ls->unconstrained) {
        PRINT1(" problem is unconstrained\n");
    }
    else {
        PRINT1("\n   [-inf,0]    [0,inf]  Finite bl  Finite bu  "
               " Two bnds      Fixed       Free\n");
        PRINT1("  %9d  %9d  %9d  %9d  %9d  %9d  %9d\n",
               nneg, npos, nlow, nupp, nnrm, nfix, nfre );
    }
    */
    return 0;

 infeas:
    //    PRINT1(" Infeasible problem: bl[%d] = %10.2e > %10.2e bu[%d]\n",
    //	   j, blj, buj, j);

    return BCLS_EXIT_INFEA;
}

/*!

  \brief Print the BCLS parameters.

  Print the user-settable BCLS parameters.
  Print level must be >= 2 for anything to print.

  \param[in] ls  BCLS solver context.

*/
void
bcls_print_params( BCLS *ls )
{
    PRINT2("\n");
    PRINT2(" Parameters:\n");
    PRINT2("  Rows............. %8d\t",   ls->m);
    PRINT2("  Columns.......... %8d\n",   ls->n);
    PRINT2("  Optimality tol... %8.1e\t", ls->optTol);
    PRINT2("  Major itn limit.. %8d\n",   ls->itnMajLim);
    PRINT2("  Damp............. %8.1e\t", ls->damp);
    PRINT2("  Minor itn limit.. %8d\n",   ls->itnMinLim);
    PRINT2("  Linsearch type... %8s\t", 
	   ls->proj_search == BCLS_PROJ_SEARCH_APPROX
	   ? "approx" : "exact");
    PRINT2("  Linear term...... %8s\n",
	   ls->c == NULL ? "no" : "yes" );
    PRINT2("  Newton step...... %8s\t",
	   ls->newton_step == BCLS_NEWTON_STEP_LSQR
	   ? "lsqr" : "cgls");
    PRINT2("  Gradient step.... %8s\n",
           ls->anorm == NULL
           ? "unscaled" : "scaled");
    PRINT2("  Preconditoner.... %8s\n",
           ls->Usolve == NULL
           ? "no" : "yes");

    return;
}

/*!

  \brief Print some stats about the column scales provided by the user.

  \param[in]  ls     BCLS solver context.
  \param[in]  anorm  Norms of the columns of A.

*/
void
bcls_examine_column_scales( BCLS *ls, double anorm[] )
{
    int j;
    int    minj;     // Index of the column with the minimum col scale.
    int    maxj;     // ............................ maximum col scale.
    double minnorm;  // The value of the minimum             col scale.
    double maxnorm;  // ........................             col scale.

    if ( anorm == NULL ) return;
    if ( ls->print_level < 2 ) return;

    minnorm = 1.0e20;
    maxnorm = 0.0;

    for (j = 0; j < ls->n; j++) {
        if ( minnorm > anorm[j] ) { minnorm = anorm[j]; minj = j; };
        if ( maxnorm < anorm[j] ) { maxnorm = anorm[j]; maxj = j; };
    }
    PRINT2( "\n" );
    PRINT2( " Column scales:  %6s   %7s\n", "column", "norm" );
    PRINT2( "       minimum   %6d   %7.1e\n", minj, minnorm  );
    PRINT2( "       maximum   %6d   %7.1e\n", maxj, maxnorm  );

    return;
}

/*!

  \brief Send output to the user-supplied pretty-printer.

  This routine prints an informative message specified by the format
  control string fmt and the optional parameter list.  See
  #bcls_set_print_hook.

  \param[in] ls   BCLS solver context.
  \param[in] current_print_level Print level of the current call.
  \param[in] fmt  String
  \param[in] ...  Any number of arguments.

*/
void
bcls_print( BCLS *ls, int current_print_level, char *fmt, ...)
{ 
    va_list arg;
    char msg[4095+1];

    // Exit immediately if the current print level is too low.
    if ( ls->print_level < current_print_level )
        return;

    // Format the message.
    va_start(arg, fmt);
    vsprintf(msg, fmt, arg);
    assert(strlen(msg) <= 4095);
    va_end(arg);
    
    // Pass the message to the user-defined hook routine.
    if (ls->print_hook != NULL && ls->print_hook(ls->print_info, msg) == 0)
        return;

    // Othwerise, send the message to the standard output */
      fprintf(stdout, "%s", msg);

      return;
}

/*!
  \brief Send output to the user-supplied error-printer and then exit.

  This routine prints an error message specified by the format control
  string fmt and the optional parameter list.  It then terminates
  execution of the program.  See #bcls_set_fault_hook.

  \param[in] ls   BCLS solver context.
  \param[in] fmt  String
  \param[in] ...  Any number of arguments.

*/
void
bcls_fault( BCLS *ls, char *fmt, ...)
{ 
    va_list arg;
    char msg[4095+1];

    // Format the message.
    va_start(arg, fmt);
    vsprintf(msg, fmt, arg);
    assert(strlen(msg) <= 4095);
    va_end(arg);
    
    // Pass the message to the user-defined hook routine.
    if (ls->fault_hook != NULL &&
	ls->fault_hook(ls->fault_info, msg) == 0) goto skip;

    // Othwerise, send the message to the standard output */
      fprintf(stderr, "%s\n", msg);

 skip:
      // return to the calling program.
      exit(EXIT_FAILURE);     
}

/*!

  \brief  Maximum primal infeasibility.

  Compute the primal infeasibility at the point x.  This is measured
  as the maximum violation of the bounds.

  \param[in]   n      Length of x, bl, bu
  \param[in]   x      Current point.
  \param[in]   bl     Lower bounds.
  \param[in]   bu     Upper bounds.
  \param[out]  pInf   Maximum infeasibility.
  \param[out]  jpInf  Variable that achieves the maximum infeasibility.

  \return
  - 0 if x is feasible to within machine precision (DBL_EPSILON).
  - 1 otherwise.

*/
int
bcls_primal_inf( int n, double x[], double bl[], double bu[],
		 double *pInf, int *jpInf )
{
    int j;
    double pj = 0;
    double blj, buj, xj;

    *pInf  = 0.0;
    *jpInf = -1;
    
    for (j = 0; j < n; j++) {
	blj = bl[j];
	buj = bu[j];
	xj  =  x[j];
	if ( xj > buj )
	    pj = xj - buj;
	else if ( blj > xj )
	    pj = blj - xj;
	
	if (*pInf  < pj ) {
	    *pInf  = pj;
	    *jpInf = j;
	}
	
    }

    return *pInf > DBL_EPSILON;
}

/*!

  \brief Maximum dual infeasibility.

  Compute the dual infeasibility (i.e., the optimality) at the point
  x, measured as the maximum departure from complementarity (with
  respect to the gradient).  If  xj  has a lower bound  lj , the dual
  infeasibility is (xj - lj ) * zj, but to allow for infinite bounds,
  we use min( xj - lj, one ) * zj.

  \param[in]  n     The length of x, z, bl, bu.
  \param[in]  x     The current point.
  \param[in]  z     The gradient at current point.
  \param[in]  bl    The lower bounds on x.
  \param[in]  bu    The upper bounds on x.
  \param[out] dInf  The maximum dual infeasibility.
  \param[out] jInf  The variable with the maximum dual infeasibility.
  
*/
void
bcls_dual_inf( int n, double x[], double z[], double bl[], double bu[],
	       double *dInf, int *jInf )
{
    int    j;
    double dj, blj, buj;

    *dInf = 0.0;
    *jInf = 0;
   
    for (j = 0; j < n; j++) {
	blj = bl[j];
	buj = bu[j];
	if (blj < buj) {
	    dj  = z[j];
	    if      (dj > 0)
		dj =   dj * fmin( x[j] - blj, 1.0 );
	    else if (dj < 0)
		dj = - dj * fmin( buj - x[j], 1.0 );
     
	    if (*dInf  < dj) {
		*dInf  = dj;
		*jInf  = j;
	    }
	}
    }

    return;
}

/*!

  \brief Project the point x into the feasible box.

  Project the point x into the box defined by the bounds bl and bu.
  
  \param[in]      n   The length of x, bl, bu.
  \param[in,out]  x   The current point.
  \param[in]      bl  The lower bounds on x.
  \param[in]      bu  The upper bounds on x.

*/
void
bcls_mid( int n, double x[], double bl[], double bu[] )
{    
    int j;

    for (j = 0; j < n; j++ ) {
	if      ( x[j] < bl[j] ) x[j] = bl[j];
	else if ( x[j] > bu[j] ) x[j] = bu[j];
    }
    return;
}

/*!

  \brief Interface to the user's matrix-vector product routine.

  This routine calls the user's Mat-vec routine, passed in as Aprod.
  We could call Aprod directly, but this interface keeps track of
  timings and call counts.

  \param[in,out]  ls    The BCLS solver context.
  \param[in]      mode  Determines if returning y = A x or x = A'y.
  \param[in]      nix   Length of ix, x, y.
  \param[in]      ix    Indices of contributing variables.
  \param[in,out]  x     RHS or LHS (depending on the mode).
  \param[in,out]  y     RHS or LHS (depending on the mode).  
  
*/
void
bcls_aprod( BCLS *ls, int mode, int nix,
	    int ix[], double x[], double y[] )
{    
    // Increase the mat-vec counters, unless this is an initialization
    // or de-initialization call.
    if (mode == BCLS_PROD_INIT ||
        mode == BCLS_PROD_TERM )
        ;// Relax.
    else {
        if      (nix == 1    ) ls->nAprod1++;
        else if (nix  < ls->n) ls->nAprodF++;
                               ls->nAprodT++;
    }

    // Call the user's routine.
    bcls_timer( &(ls->stopwatch[BCLS_TIMER_APROD]), BCLS_TIMER_START );

    int err = ls->Aprod( mode, ls->m, ls->n, nix, ix, x, y, ls->UsrWrk );

    bcls_timer( &(ls->stopwatch[BCLS_TIMER_APROD]), BCLS_TIMER_STOP );

    // Return to top of BCLS's stack if error is signaled.
    if (err)
	longjmp(ls->jmp_env, BCLS_EXIT_APROD);
    return;
}

/*!

   \brief Interface to the user's preconditioning routine.
   
   This routine calls the user's preconditioning routine, passed in as
   Usolve.  We could call Usolve directly, this this interface keeps
   tracks of timings and call counts.

   \param[in,out] ls    The BCLS solver context.
   \param[in]     mode  Determines if solving U v = w  or  U'w = v.
   \param[in]     nix   Length of ix, v, w.
   \param[in]     ix    Indices of contributing variables.
   \param[in,out] v     The solution or RHS (depending on the mode).
   \param[in,out] w     The solution or RHS (depending on the mode).

*/
void
bcls_usolve( BCLS *ls, int mode, int nix,
             int ix[], double v[], double w[] )
{    
    // Increase the counters only for actual solves.
    if (mode == BCLS_PRECON_U  ||  mode == BCLS_PRECON_Ut )
        ls->nUsolve++;

    // Call the user's routine.
    bcls_timer( &(ls->stopwatch[BCLS_TIMER_USOLVE]), BCLS_TIMER_START );

    int err = ls->Usolve( mode, ls->m, ls->n, nix, ix, v, w, ls->UsrWrk );

    bcls_timer( &(ls->stopwatch[BCLS_TIMER_USOLVE]), BCLS_TIMER_STOP );

    // Return to top of BCLS's stack if error is signaled.
    if (err)
	longjmp(ls->jmp_env, BCLS_EXIT_USOLVE);
    return;
}

/*!

   \brief Interface to the user's callback routine.
   
   This routine calls the user's callback routine.  We could call
   CallBack directly, this this interface keeps tracks of various odds
   and ends.

   \param[in,out] ls    The BCLS solver context.

   \return

   Returns (unchanged) the return from the user's callback function.

*/
int
bcls_callback( BCLS *ls )
{    
    if ( ls->CallBack )
        return ls->CallBack( ls, ls->UsrWrk );
    else
        return 0;
}

/*!

  \brief Determine the set of free variables.

  Construct an index set (stored in ix[]) of the free (i.e.,
  "inactive") variables.  The j-th variable is considered free if it
  is further than epsx from one of its boundaries, or if its gradient
  points away from that bound.

  \param[in]   ls   BCLS solver context.
  \param[in]   n    The length of x, b, bl, bu.
  \param[out]  ix   The index set of "free" variables.
  \param[in]   x    The current iterate.
  \param[in]   g    The gradient at x.
  \param[in]   bl   The lower bounds.
  \param[in]   bu   The upper bounds.

  \return

  Returns the number of free variables.  This is therefore the length
  of the constructed array ix[].

  \todo The parameter eps and epsx are used in a rather arbitrary way
  to define if g is "zero" or if x is "close" to its bound.  Need to
  understand how sensitive the algorithm is to this parameter.

*/
int
bcls_free_vars( BCLS *ls, int n, int ix[], double x[], 
		double g[], double bl[], double bu[] )
{
    int j;
    int nFree = 0;
    int low, upp, inc, dec;
    const double epsx = ls->epsx;
    const double eps  = ls->eps;

    for (j = 0; j < n; j++) {

	low  = x[j]  - bl[j]  <  epsx;
	upp  = bu[j] -  x[j]  <  epsx;

	inc  = g[j]  <  eps;
	dec  = g[j]  >  eps;

	if      ( low )
	    ; // Relax -- lower bound binding.
	else if ( upp )
	    ; // Relax -- upper bound binding.
	else {
	      // Var is free to move.
	    ix[nFree] = j;
	    nFree++;
	}
    }
    return nFree;
}

/*!
  
  \brief Discard the smallest element and contract the heap.

  On entry, the numElems of the heap are stored in x[0],...,x[numElems-1],
  and the smallest element is x[0].  The following operations are performed:
    -# Swap the first and last elements of the heap
    -# Shorten the length of the heap by one.
    -# Restore the heap property to the contracted heap.
       This effectively makes x[0] the next smallest element
       in the list.  

  \param[in]     numElems   The number of elements in the current heap.
  \param[in,out] x          The array to be modified.
  \param[in,out] ix         The indices of the array.

  \return  The number of elements in the heap after it has been contracted.

*/
int
bcls_heap_del_min( int numElems, double x[], int ix[] )
{
    assert(numElems > 0);
    
    int lastChild = numElems - 1;

    // Swap the smallest element with the lastChild.
    bcls_dswap(  x[0],  x[lastChild] );
    bcls_iswap( ix[0], ix[lastChild] );

    // Contract the heap size, thereby discarding the smallest element.
    lastChild--;
    
    // Restore the heap property of the contracted heap.
    bcls_heap_sift( 0, lastChild, x, ix );

    return numElems - 1;
}

/*!

  \brief Perform the "sift" operation for the heap-sort algorithm.

  A heap is a collection of items arranged in a binary tree.  Each
  child node is less than or equal to its parent.  If x[k] is the
  parent, than its children are x[2k+1] and x[2k+2].

  This routine promotes ("sifts up") children that are smaller than
  their parents.  Thus, this is a "reverse" heap, where the smallest
  element of the heap is the root node.

  \param[in]     root       The root index from which to start sifting.
  \param[in]     lastChild  The last child (largest node index) in the sift operation.
  \param[in,out] x          The array to be sifted.
  \param[in,out] ix         The indices of the array.

*/
void
bcls_heap_sift( int root, int lastChild, double x[], int ix[] )
{
    int child;

    for (; (child = (root * 2) + 1) <= lastChild; root = child) {

	if (child < lastChild)
	    if ( x[child] > x[child+1] )
		child++;
	
	if ( x[child] >= x[root] )
	    break;

	bcls_iswap( ix[root], ix[child] );
	bcls_dswap(  x[root],  x[child] );
    }
    return;
}

/*!
  
  \brief  Build a heap by adding one element at a time.
  
  \param[in]      n   The length of x and ix.
  \param[in,out]  x   The array to be heapified.
  \param[in,out]  ix  The indices of the array elements.

*/
void
bcls_heap_build( int n, double x[], int ix[] )
{    
    int i;

    for (i = n/2; i >= 0; i--)
	bcls_heap_sift( i, n-1, x, ix );

    return;
}

/*!

  \brief Return the 2-norm of x(ix), where ix is a set of indices.
  
  \param[in] nix  Length of ix.
  \param[in] ix   Set of indices to access in x.
  \param[in] x    Input vector.
  
  \return         The 2-norm of the subvector.

*/
double
bcls_vec_dnorm2( int nix, int ix[], double x[] )
{    
    int j, k;
    double norm = 0.0;
    
    for (j = 0; j < nix; j++) {
        k     = ix[j];
        norm +=  x[k] * x[k];
    }
    norm = sqrt(norm);
    
    return norm;
}

/*!

  \brief Load a constant into a vector.

  Loads a constant into every component of a vector x.
  The (usual) case incx = 1, alpha = 0.0 is treated specially.

  \param[in]      n      The length of x.
  \param[in]      alpha  The constant double to be loaded into x.
  \param[in,out]  x      The vector to be loaded.
  \param[in]      incx   Fill every "incx" elements.

*/
void
bcls_dload( const int n, const double alpha, double x[], const int incx )
{    
    int i, ix;

    if (incx <= 0) return;

    ix = 0;

    for (i = 0; i < n; i++) {
        x[ix]  = alpha;
        ix    += incx;
    }
    return;
}

/*!

  \brief Gather the components of x(ix) into y, eg, y = x(ix).

  \param[in]     nix   The length of the index set.
  \param[in]     ix    The index set.
  \param[in]     x     The vector to be gathered.
  \param[out]    y     The gathered vector.

*/
void
bcls_gather( const int nix, const int ix[], const double x[], double y[] )
{
    int j, k;
    for (j = 0; j < nix; j++) {
        k = ix[ j ];
        y[ j ] = x[ k ];
    }
}

/*!

  \brief Scatter the components of x into y(ix), eg, y(ix) = x.

  \param[in]     nix   The length of the index set.
  \param[in]     ix    The index set.
  \param[in]     x     The vector to be scattered.
  \param[out]    y     The scattered vector.

*/
void
bcls_scatter( const int nix, const int ix[], const double x[], double y[] )
{
    int j, k;
    for (j = 0; j < nix; j++) {
        k = ix[ j ];
        y[ k ] = x[ j ];
    }
}
