#
# Makefile for fsgdf
#

if ENABLE_TCLTK_APPS
AM_CPPFLAGS=-I$(top_srcdir)/include $(TCL_CFLAGS)
else
AM_CPPFLAGS=-I$(top_srcdir)/include
endif
AM_LDFLAGS=

noinst_LIBRARIES= libfsgdf.a
libfsgdf_a_SOURCES = fsgdf.c

if ENABLE_TCLTK_APPS
check_PROGRAMS=test test_wrapper
test_SOURCES=test.c fsgdf_wrap.c
test_LDFLAGS=$(LDFLAGS_TCL)
test_LDADD=libfsgdf.a $(addprefix $(top_srcdir)/, $(LIBS_MGH)) $(LIBS_TCL)
TESTS=test_wrapper
test_wrapper:
	cp $(top_srcdir)/scripts/tcltktixblt_wrapper test_wrapper
	echo "setenv SURFER_FRONTDOOR" >> test_wrapper
	echo "./test \$$argv" >> test_wrapper
	chmod a+x test_wrapper
endif

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
