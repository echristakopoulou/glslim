include $(top_srcdir)/Makefile.common

AUTOMAKE_OPTIONS = foreign

INCLUDES = -I$(top_srcdir)/include

lib_LTLIBRARIES = \
	libbcls.la \
	libcblas.la

libbcls_la_SOURCES = \
	bcls.c \
	bclib.c \
	bclsqr.c \
	bccgls.c \
	bcsolver.c \
	bctimer.c \
	bcversion.c \
	lsqr.c
if ENABLE_DEBUGGING
libbcls_la_SOURCES += bcdebug.c
endif

libcblas_la_SOURCES = \
	bccblas.c

include_HEADERS = \
	$(top_srcdir)/include/bcls.h

EXTRA_DIST = \
	$(top_srcdir)/include/bccgls.h \
	$(top_srcdir)/include/bclib.h \
	$(top_srcdir)/include/bclsqr.h \
	$(top_srcdir)/include/bcsolver.h \
	$(top_srcdir)/include/bctimer.h \
	$(top_srcdir)/include/bcversion.h \
	$(top_srcdir)/include/cblas.h \
	$(top_srcdir)/include/lsqr.h

libcblas_la_LDFLAGS = -version-info $(BCLS_LIBRARY_VERSION)
libbcls_la_LDFLAGS  = -version-info $(BCLS_LIBRARY_VERSION)

debug:
	$(MAKE) CFLAGS="$(CFLAGS_DEBUG)"


