include $(top_srcdir)/Makefile.common

AM_CPPFLAGS = -I$(top_srcdir)/include

# Less strict checking.
AUTOMAKE_OPTIONS = foreign

noinst_PROGRAMS = bcsol
bcsol_SOURCES   = bcsol.c cs_util.c cs_malloc.c iohb.c iohb.h cs.h
bcsol_LDADD     = $(top_builddir)/src/libbcls.la
if ENABLE_LOCAL_CBLAS
bcsol_LDADD    += $(top_builddir)/src/libcblas.la
else
bcsol_LDADD    += $(CBLAS_LIBS) $(BLAS_LIBS)
endif

# Other things to distribute.
EXTRA_DIST = \
	mat2hb.m \
	mat2hbmex.c \
	savehb.m \
	dm2hb.m \
	illc1033.rra \
	illc1850.rra \
	well1033.bnd \
	well1033.rra \
	well1850.rra

debug:
	$(MAKE) CFLAGS="$(CFLAGS_DEBUG)"

if HAVE_MATLAB

AM_CPPFLAGS += -I$(top_srcdir) -I$(MATLAB_DIR)/extern/include

mat2hbmex.o: mat2hbmex.c
	$(LTCOMPILE) -c -o $@ $<

mex: mat2hbmex.$(MEX_EXT)

mat2hbmex.$(MEX_EXT): mat2hbmex.o iohb.o
	$(LIBTOOL) --tag=CC --mode=link $(MEX) -static -no-install -cxx -o mat2hbmex $^
endif

clean-local:
	-rm -f mat2hbmex.mex*
	-rm -f *.o
