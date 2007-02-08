##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mris_compute_optimal_kernel
mris_compute_optimal_kernel_SOURCES=mris_compute_optimal_kernel.c
mris_compute_optimal_kernel_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_compute_optimal_kernel_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
