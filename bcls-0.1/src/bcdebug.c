/* bcdebug.c
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
   Debugging routines.  Only linked into library during "make debug".
*/

#include "bclib.h"

unsigned long long int hashd(double x)
{
    union 
    { 
	unsigned long long int i;
	double d; 
    } u;
    if (x == 0.0)
        return 0ULL;
    u.d = x;
    return u.i;
}

unsigned long long int hashiv(const int n, int x[])
{
    unsigned long long int sum = 0ULL;
    int j;
    if (!x)
        return sum;
    for (j = 0; j < n; j++) {
        if (x[j] == 0)
            ; // Relax.
        else
            sum += (unsigned long long int)x[j];
    }
    return sum;
}

unsigned long long int hashdv(const int n, double x[])
{
    union 
    { 
	unsigned long long int i;
	double d; 
    } u;
    unsigned long long int sum = 0ULL;
    int j;
    if (!x)
        return sum;
    for (j = 0; j < n; j++) {
        if (x[j] == 0.0)
            ; // Relax.
        else {
            u.d  = x[j];
            sum += u.i;
        }
    }
    return sum;
}

unsigned long long int hashls( BCLS *ls )
{
    unsigned long long int sum = 0ULL;
    int
        n  = ls->n,
        m  = ls->m,
        nm = n + m;

    sum += n + m + nm;
    sum += hashdv(n  , ls->anorm);
    sum += hashdv(n  , ls->x);
    sum += hashdv(m  , ls->b);
    sum += hashdv(n  , ls->c);
    sum += hashdv(n  , ls->bl);
    sum += hashdv(n  , ls->bu);
    sum += hashdv(nm , ls->r);
    sum += hashdv(n  , ls->g);
    sum += hashdv(n  , ls->dx);
    sum += hashdv(n  , ls->dxFree);
    sum += hashdv(n  , ls->aBreak);
    sum += hashiv(n  , ls->iBreak);
    sum += hashiv(n  , ls->ix);
    sum += hashdv(nm , ls->wrk_u);
    sum += hashdv(nm , ls->wrk_v);
    sum += hashdv(nm , ls->wrk_w);

    return sum;
}
