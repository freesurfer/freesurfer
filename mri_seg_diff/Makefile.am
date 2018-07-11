##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mri_seg_diff
mri_seg_diff_SOURCES=mri_seg_diff.c
mri_seg_diff_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_seg_diff_LDFLAGS=$(OS_LDFLAGS)

TESTS=test.py

EXTRA_DIST= test.py

EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
