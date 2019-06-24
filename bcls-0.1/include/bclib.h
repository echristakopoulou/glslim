/* bclib.h
   $Revision: 282 $ $Date: 2006-12-17 17:38:00 -0800 (Sun, 17 Dec 2006) $

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

#ifndef _BCLSLIB_H
#define _BCLSLIB_H

#include "bcls.h"

#define PRINT1(...)  bcls_print( ls, 1, __VA_ARGS__ )
#define PRINT2(...)  bcls_print( ls, 2, __VA_ARGS__ )
#define PRINT3(...)  bcls_print( ls, 3, __VA_ARGS__ )
#define PRINT4(...)  bcls_print( ls, 4, __VA_ARGS__ )
#define PRINT5(...)  bcls_print( ls, 5, __VA_ARGS__ )
#define PRINT6(...)  bcls_print( ls, 6, __VA_ARGS__ )

// Inline min and max functions.
static inline int imax(int a, int b) {
    return a > b ? a : b;
}

static inline int imin(int a, int b) {
    return a > b ? b : a;
}

// Inline int swap function.
#define bcls_iswap(a,b) {	\
        int c = (a);    \
        (a) = (b);      \
        (b) =  c;       \
    }

// Inline double swap function.
#define bcls_dswap(a,b) {    \
	double c = (a);	\
	(a) = (b);      \
	(b) =  c;       \
    }

// Print output.
void bcls_print( BCLS *ls, int current_print_level, char *fmt, ... );

// Print fault and exit the program.
void bcls_fault( BCLS *ls, char *fmt, ... );

// Compute a projected step.
void bcls_project_step( int n, double s[], double step, double x[],
                        double dx[], double bl[], double bu[] );

// Examine and print something about the problem's bounds.
int bcls_examine_bnds( BCLS *ls, int n, double bl[], double bu[] );

// Examine and print something about the column scales.
void bcls_examine_column_scales( BCLS *ls, double anorm[] );

// Print the BCLS parameters.
void bcls_print_params( BCLS *ls );

// Compute the primal infeasibility at the point x.
int bcls_primal_inf( int n, double x[], double bl[], double bu[],
                     double *pInf, int *jpInf );

// Compute the dual infeasibility at the point x.
void bcls_dual_inf( int n, double x[], double z[], double bl[], double bu[],
		    double *dInf, int *jInf );

// Project  x  into the box defined by bl and bu.
void bcls_mid( int n, double x[], double bl[], double bu[] );

// Interface to the user's Aprod routine.
void bcls_aprod( BCLS *ls, int mode, int nix,
		 int ix[], double x[], double y[] );

// Interface to the user's preconditioning routine.
void bcls_usolve( BCLS *ls, int mode, int nix,
                  int ix[], double v[], double w[] );

// Interface to the user's callback routine.
int bcls_callback( BCLS *ls );

// Determine the set of free variables.
int bcls_free_vars( BCLS *ls, int n, int ix[], double x[],
		    double g[], double bl[], double bu[] );

// Discard the smallest element from the heap.
int bcls_heap_del_min( int numElems, double x[], int ix[] );

// Sift a heap.
void bcls_heap_sift( int root, int lastChild, double x[], int ix[] );

// Build a heap.
void bcls_heap_build( int n, double x[], int ix[] );

// Two-norm of a subvector.
double
bcls_vec_dnorm2( int nix, int ix[], double x[] );

// Load a constant into a vector.
void bcls_dload( const int n, const double alpha, double x[], const int incx );

// Gather the components of x(ix) into y, eg, y = x(ix).
void bcls_gather( const int nix, const int ix[],
                  const double x[], double y[] );

// Scatter the components of x into y(ix), eg, y(ix) = x.
void bcls_scatter( const int nix, const int ix[],
                   const double x[], double y[] );

#endif
