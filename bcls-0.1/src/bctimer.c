/* bctimer.c
   $Revision: 237 $ $Date: 2006-04-19 18:42:26 -0700 (Wed, 19 Apr 2006) $

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
   Timer functions.
*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <assert.h>
#include "bctimer.h"

/*!
  \brief Start and stop the various stopwatches.
  
  \param[in,out] swatch  The stopwatch operated on.
  \param[in]     task    Start, stop, init, or print the stopwatch.

  \return  Total time on the stopwatch.

*/
double
bcls_timer( BCLS_timer *swatch, int task )
{
    assert( task == BCLS_TIMER_START ||
	    task == BCLS_TIMER_STOP  ||
	    task == BCLS_TIMER_INIT  ||
	    task == BCLS_TIMER_PRINT );

    double last; // Wall-clock time that the stopwatch was last called.

#ifdef HAVE_GETRUSAGE

    struct rusage ru;
    (void) getrusage(RUSAGE_SELF, &ru);

    double user_time = 
                 ru.ru_utime.tv_sec   // user time (seconds)
        + 1e-6 * ru.ru_utime.tv_usec; // user time (microseconds)

    double sys_time = 
                 ru.ru_stime.tv_sec   // system time (seconds)
        + 1e-6 * ru.ru_stime.tv_usec; // system time (microseconds)

    last = user_time + sys_time;

#else

    last = (double)clock();

#endif

    if ( task == BCLS_TIMER_START ) {
	swatch->start = last;
	swatch->nStarts++;
    }
    else if ( task == BCLS_TIMER_STOP ) {

#ifdef HAVE_GETRUSAGE

	swatch->total += last - swatch->start;

#else

        swatch->total += (last - swatch->start)
                             / ((double)CLOCKS_PER_SEC);

#endif

    }
    else if ( task == BCLS_TIMER_INIT ) {
	swatch->total   = 0.0;
	swatch->nStarts = 0;
    }
    else if ( task == BCLS_TIMER_PRINT ) {
	return swatch->total;
    }

    return 0.0;
}
