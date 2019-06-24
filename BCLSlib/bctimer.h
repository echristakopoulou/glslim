/* bctimer.h
   $Revision: 238 $ $Date: 2006-04-19 18:42:33 -0700 (Wed, 19 Apr 2006) $

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
   BCLS timer functions.
*/

#ifndef _BCLSTIMER_H
#define _BCLSTIMER_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_GETRUSAGE
#include <sys/time.h>
#include <sys/resource.h>
#else
#include <time.h>
#endif

/* BCLS Timer */
typedef struct BCLS_timer BCLS_timer;
struct BCLS_timer {
    double  start;
    double  total;
    int     nStarts;
    char   *name;
};

/* Measure elapsed time between calls. */
#define BCLS_TIMER_INIT  -1
#define BCLS_TIMER_START  0
#define BCLS_TIMER_STOP   1
#define BCLS_TIMER_PRINT  2
double bcls_timer( BCLS_timer *timer, int task );

#endif
