include $(top_srcdir)/Makefile.common

# Less strict checking
AUTOMAKE_OPTIONS = no-dependencies foreign

AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/include -I$(MATLAB_DIR)/extern/include

# Other things to distribute.
EXTRA_DIST = \
	aprod.m \
	bcls.m \
	bclsmex.c \
	cblas.c \
	test.m \
	usolve.m \
	well1033.mat

# Create a library  This convinces Automake to setup the
# LIBTOOL and LTCOMPILE variables used below.
noinst_PROGRAMS = foo
foo_SOURCES = foo.c

mex: bclsmex.$(MEX_EXT)

if HAVE_MATLAB

.c.lo:
	$(LTCOMPILE) -c -o $@ $<

bclsmex.$(MEX_EXT): bclsmex.lo cblas.lo
	$(LIBTOOL) --tag=CC --mode=link $(MEX) -static -no-install -cxx -o bclsmex $(top_builddir)/src/libbcls.la $^

endif

clean-local:
	-rm -f bclsmex.mex*
	-rm -f *.o
