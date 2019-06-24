/* cblas.c
   $Revision$ $Date$

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

   This file contains C-wrappers to the BLAS (Basic Linear Algebra
   Subprograms) routines that are used by BCLS.  It's intended only
   to be used with compiling the MEX interface to BCLS, and makes
   use of the BLAS routines supplied by MATLAB.  Each routine is a
   direct call to the Fortran verision of that BLAS subroutine.

   Included BLAS routines:

   - cblas_daxpy
   - cblas_dcopy
   - cblas_ddot
   - cblas_dnrm2
   - cblas_dscal
   - cblas_dgemv
*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "cblas.h"

#ifdef WINDOWS_HOST
#define F77_FUNC(name,NAME) name
#else
#define F77_FUNC(name,NAME) name ## _
#endif

/*!
  \param[in]     N
  \param[in]     alpha
  \param[in]     X      
  \param[in]     incX
  \param[in,out] Y
  \param[in]     incY
*/
void
cblas_daxpy( const int N, const double alpha, const double *X,
             const int incX, double *Y, const int incY)
{
    void F77_FUNC(daxpy,DAXPY)();
    F77_FUNC(daxpy,DAXPY)(&N, &alpha, X, &incX, Y, &incY);
    return;
}

/*!
  \param[in]     N
  \param[in]     X      
  \param[in]     incX
  \param[out]    Y
  \param[in]     incY
*/
void
cblas_dcopy( const int N, const double *X,
             const int incX, double *Y, const int incY)
{
    void F77_FUNC(dcopy,DCOPY)();
    F77_FUNC(dcopy,DCOPY)(&N, X, &incX, Y, &incY);
    return;
}

/*!
  \param[in]     N
  \param[in]     X      
  \param[in]     incX
  \param[in]     Y
  \param[in]     incY
  
  \return  Dot product of X and Y.

*/
double
cblas_ddot( const int N, const double *X,
            const int incX, const double *Y, const int incY)
{
    double F77_FUNC(ddot,DDOT)();
    return F77_FUNC(ddot,DDOT)(&N, X, &incX, Y, &incY);
}

/*!
  \param[in]     N
  \param[in]     X      
  \param[in]     incX

  \return Two-norm of X.
*/
double
cblas_dnrm2( const int N, const double *X, const int incX) 
{
    double F77_FUNC(dnrm2,DNRM2)();
    return F77_FUNC(dnrm2,DNRM2)(&N, X, &incX);
}

/*!
  \param[in]     N
  \param[in]     alpha
  \param[in,out] X
  \param[in]     incX
*/
void
cblas_dscal(const int N, const double alpha, double *X, const int incX)
{
    void F77_FUNC(dscal,DSCAL)();
    F77_FUNC(dscal,DSCAL)(&N, &alpha, X, &incX);
    return;
}


void
cblas_dgemv(const enum CBLAS_ORDER order,
            const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
            const double alpha, const double  *A, const int lda,
            const double  *X, const int incX, const double beta,
            double  *Y, const int incY)
{
    void F77_FUNC(dgemv,DGEMV)();
    char TA;
    if (order == CblasColMajor) {
        if      (TransA == CblasNoTrans)   TA = 'N';
        else if (TransA == CblasTrans)     TA = 'T';
        else if (TransA == CblasConjTrans) TA = 'C';
        F77_FUNC(dgemv,DGEMV)
            (&TA, &M, &N, &alpha, A, &lda, X, &incX, &beta, Y, &incY);
    }
    else if (order == CblasRowMajor) {
        if      (TransA == CblasNoTrans)   TA = 'T';
        else if (TransA == CblasTrans)     TA = 'N';
        else if (TransA == CblasConjTrans) TA = 'N';
        F77_FUNC(dgemv,DGEMV)
            (&TA, &N, &M, &alpha, A, &lda, X, &incX, &beta, Y, &incY);
    }
    return;
}

