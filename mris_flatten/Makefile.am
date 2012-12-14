## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

bin_PROGRAMS = mris_flatten
mris_flatten_SOURCES=mris_flatten.c
mris_flatten_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_flatten_LDFLAGS=$(OS_LDFLAGS)

## ----
## CUDA
## ----

# BUILDCUDA is defined if configure.in finds CUDA
if BUILDCUDA

NVCCFLAGS += -g -arch sm_13 --ptxas-options=-v

# rules for building cuda files
.cu.o:
	$(NVCC) -o $@ -c $< $(NVCCFLAGS) $(CUDA_CFLAGS) -I$(top_srcdir)/include

bin_PROGRAMS += mris_flatten_cuda
mris_flatten_cuda_SOURCES = mris_flatten.c
mris_flatten_cuda_CFLAGS = $(CUDA_CFLAGS) $(AM_CFLAGS) -DFS_CUDA
mris_flatten_cuda_CXXFLAGS = $(CUDA_CFLAGS) $(AM_CFLAGS) -DFS_CUDA
mris_flatten_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS)
mris_flatten_cuda_LDFLAGS = $(OS_LDFLAGS) 
mris_flatten_cuda_LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(mris_flatten_cuda_CFLAGS) \
        $(CFLAGS) $(mris_flatten_cuda_LDFLAGS) $(LDFLAGS) -o $@

endif

EXTRA_DIST = make_mris_flatten_movie.tcl

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
