/* bcversion.c
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
   Record the version and latest compliation of the BCLS library.
*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#define PACKAGE_VERSION "NN"
#endif
#include "bcversion.h"

static char *version = PACKAGE_VERSION;
static char *last_compilation = __DATE__ " " __TIME__;

/*!
  Return the version of the library.
*/
char *
bcls_version_info( )
{
    return version;
}

/*!
  Return a timestamp of the last compilation.
*/
char *
bcls_compilation_info( )
{
    return last_compilation;
}
