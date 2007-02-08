#
# Makefile.am for libvtkutils.a
#

AM_CPPFLAGS= -Wno-deprecated \
	$(X_CFLAGS) \
	$(TCL_CFLAGS) \
	$(VXL_CFLAGS) \
	-I$(top_srcdir)/include \
	$(VTK_CXXFLAGS) \
	$(KWWIDGETS_CXXFLAGS)

AM_LDFLAGS=

noinst_LIBRARIES= libvtkutils.a

libvtkutils_a_SOURCES=$(wildcard *.cxx)

am_libvtkutils_a_OBJECTS=$(libvtkutils_a_SOURCES:%.cxx=%.o)

.cxx.o:
	@echo "compiler: $(COMPILE) end_compiler"
	source='$<' object='$@' libtool=no \
	depfile='$(DEPDIR)/$*.Po' tmpdepfile='$(DEPDIR)/$*.TPo' \
	$(CCDEPMODE) $(depcomp) \
	$(COMPILE) -c `test -f '$<' || echo '$(srcdir)/'`$<

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
