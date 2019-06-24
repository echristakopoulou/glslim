/* bcls.h
   $Revision: 276 $ $Date: 2006-12-09 21:00:00 -0800 (Sat, 09 Dec 2006) $

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
   BCLS library header file.
   This is the file "included" by the user.
*/

#ifndef _BCLS_H
#define _BCLS_H

/* Prevent C++ programs from name mangling these definitions. */
#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <setjmp.h>
#include "bctimer.h"

/* The BCLS objects define everything about a problem. */
typedef struct BCLS BCLS;
struct BCLS {

    /*!\brief Transit pointer passed (untouched) to the routine print_hook. */
    void *print_info;

    /*!\brief User-defined print-hook routine. */
    int (*print_hook)(void *info, char *msg);

    /*!\brief Transit point passed (untouched) to the routine fault_hook. */
    void *fault_info;

    /*!\brief User-defined fault-hook routine. */
    int (*fault_hook)(void *info, char *msg);

    /*!\brief User-defined matrix-vector multiply routine. */
    int (*Aprod)( int mode, int m, int n, int nix,
		  int ix[], double x[], double y[], void *UsrWrk );
#define BCLS_PROD_A     1 /*!< Mode for Aprod: y = A x. */
#define BCLS_PROD_At    2 /*!< Mode for Aprod: x = A'y. */
#define BCLS_PROD_INIT -1 /*!< Mode for Aprod: Initialize. */
#define BCLS_PROD_TERM -2 /*!< Mode for Aprod: Destroy. */

    /*!\brief User-defined preconditioning routine. */
    int (*Usolve)( int mode, int m, int n, int nix,
                   int ix[], double v[], double w[], void *UsrWrk );
#define BCLS_PRECON_U     1 /*!< Mode for Usolve: U v = w. */
#define BCLS_PRECON_Ut    2 /*!< Mode for Usolve: U'w = v. */
#define BCLS_PRECON_INIT -1 /*!< Mode for Usolve: Initialize. */
#define BCLS_PRECON_TERM -2 /*!< Mode for Usolve: Destroy. */

    /*!\brief User-defined callback function. */
    int (*CallBack)( BCLS *ls, void *UsrWrk );

    /*!\brief User-defined column weights.  Length n.  Optional. */
    double *anorm;
#define BCLS_MIN_COLUMN_NORM  1e-8 /*!<\brief Min allowable col norm.*/

    /*!\brief User workspace.  Passed (untouched) to Aprod. */
    void *UsrWrk;

    /* User options. */
    int print_level;               /*!<\brief Log verbosity. */
    int proj_search;               /*!<\brief Type of linesearch. */
#define BCLS_PROJ_SEARCH_EXACT   0 /*!<\brief Exact projected linesearch. */
#define BCLS_PROJ_SEARCH_APPROX  1 /*!<\brief Approximate projected linesearch. */
    int newton_step;               /*!<\brief Method for computing Newton step. */
#define BCLS_NEWTON_STEP_LSQR    0 /*!<\brief LSQR on least-squares formulation. */
#define BCLS_NEWTON_STEP_CGLS    1 /*!<\brief CGLS on normal equations. */
    FILE *minor_file;              /*!<\brief Output file for minor iteration log.*/

    /*!<\brief Timers. */
    BCLS_timer stopwatch[4];
#define BCLS_TIMER_TOTAL   0 /*!<\brief Overall tiemr. */
#define BCLS_TIMER_APROD   1 /*!<\brief Mat-vec routine timer. */
#define BCLS_TIMER_LSQR    2 /*!<\brief LSQR timer. */
#define BCLS_TIMER_USOLVE  3 /*!<\brief Usolve timer. */

    /* Counters (iterations, no. of mat-vecs calls, etc.) */
    int itnMaj;     /*!<\brief No. of major iterations. */
    int itnMajLim;  /*!<\brief Max no. of major iterations. */
    int itnMin;     /*!<\brief No. of cummulative minor (LSQR) iterations. */
    int itnMinLim;  /*!<\brief Max no. of minor (LSQR) iterations. */
    int nAprodT;    /*!<\brief No. of mat-vec calls, with     len(ix) = n. */
    int nAprodF;    /*!<\brief No. of mat-vec calls, with 1 < len(ix) < n. */
    int nAprod1;    /*!<\brief No. of mat-vec calls, with 1 = len(ix).     */
    int nUsolve;    /*!<\brief No. of call to Usolve. */

    /*!\brief Number of rows in A. */
    int m;

    /*!\brief Number of columns in A. */
    int n;

    /*!\brief Maximum no. of rows in A allocated for. */
    int mmax;

    /*!\brief Maximum no. of columns in A allocated for. */
    int nmax;

    /*!\brief Properties of the problem. 0 if constrained; 1 if unconstrained.*/
    int unconstrained;

    /* Damping parameters. */
    double damp;            /*!<\brief Set by user.  Must be >= 0. */
    double damp_min;        /*!<\brief Only used when there's a linear term. */
    double damp_actual;     /*!<\brief Set in newton_step, used in aprod_free. */

    /*!\brief Exit condition. */
    int exit;
#define BCLS_EXIT_CNVGD    0 /*!<\brief Optimal solution found. */
#define BCLS_EXIT_MAJOR    1 /*!<\brief Too many major iterations. */
#define BCLS_EXIT_MINOR    2 /*!<\brief Too many minor iterations. */
#define BCLS_EXIT_UNDEF    3 /*!<\brief Exit condition is undefined. */
#define BCLS_EXIT_UNBND    4 /*!<\brief Found direction of infinite descent. */
#define BCLS_EXIT_INFEA    5 /*!<\brief Bounds are inconsistent. */
#define BCLS_EXIT_LFAIL    6 /*!<\brief Linesearch failure. */
#define BCLS_EXIT_APROD  100 /*!<\brief Exit requested by user's aprod  routine. */
#define BCLS_EXIT_USOLVE 110 /*!<\brief Exit requested by user's usolve routine. */
#define BCLS_EXIT_CALLBK 120 /*!<\brief Exit requested by user's callback routine. */
    
    /* Status of the solution (stored in x). */
    double soln_obj;       /*!<\brief Objective value. */
    double soln_rNorm;     /*!<\brief Residual at solution. */
    double soln_dInf;      /*!<\brief Dual infeasibility of solution. */
    int    soln_jInf;      /*!<\brief Variable with maximum dual infeasibility. */
    int    soln_stat;      /*!<\brief Status of the solution. */
#define BCLS_SOLN_UNDEF  0 /*!<\brief Solution is undefined. */
#define BCLS_SOLN_OPTIM  1 /*!<\brief Solution is optimal. */
    
    /* Tolerances and limits. */
    double optTol;          /*!<\brief Optimality tolerance. */
    double conlim;          /*!<\brief Upper limit on cond. no. of [A; gamma I]. */
    double mu;              /*!<\brief Armijo sufficient decrease factor. */
    double backtrack;       /*!<\brief Reduction factor for Armijo linesearch. */
    int    backtrack_limit; /*!<\brief Limit on no. of backtracks in linesearch. */

    /* Constants. */
    double eps;             /*!<\brief Machine precision. */
    double eps2;            /*!<\brief eps ^ 0.50. */
    double eps3;            /*!<\brief eps ^ 0.75. */
    double epsx;            /*!<\brief Tolerance for x "on bound". */
    double epsfixed;        /*!<\brief Tolerance for x considered "fixed". */
    double BigNum;          /*!<\brief Biggest number we'll allow in x. */
#define BCLS_INFINITY 1e+20

    /*!\brief x(n) solution (allocated by user). */
    double *x;

    /*!\brief b(m): right-hand-side vector (allocated by user). */
    double *b;

    /*!\brief c(n): linear term (allocated by user). */
    double *c;

    /*!\brief bl(n): lower bounds (allocated by user). */
    double *bl;

    /*!\brief bu(n): upper bounds (allocated by user). */
    double *bu;

    /* BCLS workspace (allocated by BCLS). */
    double *r;              /*!<\brief (max(n, m)) Residual. */
    double *g;              /*!<\brief (n) Gradient. */
    double *dx;             /*!<\brief (n) Search direction, full space. */
    double *dxFree;         /*!<\brief (n) Search direction, subspace. */
    double *aBreak;         /*!<\brief (n) Step to each breakpoint. */
    int    *iBreak;         /*!<\brief (n) Indices of each breakpoint. */
    int    *ix;             /*!<\brief (n) Indices of free variables. */
    double *wrk_u;          /*!<\brief (max(n, m)) Workspace. */
    double *wrk_v;          /*!<\brief (max(n, m)) Workspace. */
    double *wrk_w;          /*!<\brief (max(n, m)) Workspace. */

    /*!\brief Long jump environment. */
    jmp_buf jmp_env;

};

/* Create a new BCLS problem. */
BCLS *
bcls_create_prob( int mmax, int nmax );

/* Initialize  a new BCLS problem. */
void
bcls_init_prob( BCLS *ls );

/* Free a BCLS problem. */
int
bcls_free_prob( BCLS *ls );

/* Return an exit message. */
char *
bcls_exit_msg( int flag );

/* Driver for the BCLS solver. */
int
bcls_solve_prob( BCLS *ls );

/* Install a user-defined print-hook routine. */
void
bcls_set_print_hook( BCLS *ls, void *info,
                     int (*hook)(void *info, char *msg) );

/* Install a user-defined print-hook routine. */
void
bcls_set_fault_hook( BCLS *ls, void *info,
                     int (*hook)(void *info, char *msg) );

/* Install a user-defined preconditioning routine. */
void
bcls_set_usolve( BCLS *ls,
                 int (*Usolve)( int mode, int m, int n, int nix,
                                int ix[], double v[], double w[],
                                void *UsrWrk ) );

/* Give BCLS access to a set of column norms of A. */
void
bcls_set_anorm( BCLS *ls, double anorm[] );

/* Compute column norms of A. */
int
bcls_compute_anorm( BCLS *ls, int n, int m,
                    int (*Aprod)
                    ( int mode, int m, int n, int nix,
                      int ix[], double x[], double y[], void *UsrWrk ),
                    void *UsrWrk,
                    double anorm[] );

/* Set the problem data. */
void
bcls_set_problem_data( 
     BCLS  *ls,    /* A BCLS problem */
     int    m,     /* No. of problem rows */
     int    n,     /* No. of problem columns */
     int (*Aprod)  /* Matrix-vector product routine */
         ( int mode, int m, int n, int nix,
	   int ix[], double x[], double y[], void *UsrWrk ),
     void  *UsrWrk,/* Arbitrary user data passed to Aprod */
     double damp,  /* The damping parameter */
     double x[],   /* The solution vector (n) */
     double b[],   /* The RHS vector (m) */
     double c[],   /* The linear term (n) */
     double bl[],  /* The lower bounds vector */
     double bu[]   /* The upper bounds vector */
     );
			 
#ifdef __cplusplus
}
#endif

#endif /* _BCLS_H */
