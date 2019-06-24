/* bcsolver.h
   $Revision: 224 $ $Date: 2006-04-12 16:19:59 -0700 (Wed, 12 Apr 2006) $

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
   BCLS solver library.  All of this is private to the BCLS library.
*/

#ifndef _BCLSSOLVER_H
#define _BCLSSOLVER_H

#include "bcls.h"

// The actual BCLS solver.
void
bcls_solver( BCLS *ls, int m, int n, double *bNorm,
	     double x[], double b[], double c[],
	     double bl[], double bu[],
	     double r[], double g[], double dx[], double dxFree[],
	     int ix[], double aBreak[], int iBreak[], double anorm[] );

#endif
