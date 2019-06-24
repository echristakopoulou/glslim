include $(top_srcdir)/Makefile.common

# Less strict checking
AUTOMAKE_OPTIONS = foreign

# Directories to build.
SUBDIRS = src examples doc
if HAVE_MATLAB
SUBDIRS += matlab
endif

# By default, also make the MEX interface if Matlab was found.
if HAVE_MATLAB
all-local: mex
endif

html-doc:
	$(MAKE) -C doc html-doc

debug:
	$(MAKE) CFLAGS="$(CFLAGS_DEBUG)"

local-install:
	$(MAKE) install prefix="@abs_top_builddir@/$(build)"

mex:
	$(MAKE) $(AM_MAKEFLAGS) -C examples mex
	$(MAKE) $(AM_MAKEFLAGS) -C matlab   mex

honk:
	@echo "Honk!"