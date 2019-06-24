include $(top_srcdir)/Makefile.common

# Less strict checking.
AUTOMAKE_OPTIONS = foreign

# Other things to distribute.
EXTRA_DIST = \
	Doxyfile.in \
	license-stub.txt \
	versionfilter.sh \
	local/bcls_footer.html \
	local/bcls_header.html \
	local/index.doc \
	local/installation.doc \
	local/autoconf_install.doc


html-doc:
	doxygen
