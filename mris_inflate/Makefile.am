## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS = mris_inflate
mris_inflate_SOURCES=mris_inflate.c
mris_inflate_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_inflate_LDFLAGS=$(OS_LDFLAGS)

## ----
## CUDA
## ----

# BUILDCUDA is defined if configure.in finds CUDA
if BUILDCUDA

NVCCFLAGS += -g -arch sm_13 --ptxas-options=-v

# rules for building cuda files
.cu.o:
	$(NVCC) -o $@ -c $< $(NVCCFLAGS) -I$(top_srcdir)/include

bin_PROGRAMS += mris_inflate_cuda
mris_inflate_cuda_SOURCES = mris_inflate.c
mris_inflate_cuda_CFLAGS = $(AM_CFLAGS)  $(CUDA_CFLAGS) -DFS_CUDA
mris_inflate_cuda_CXXFLAGS = $(AM_CFLAGS) $(CUDA_CFLAGS) -DFS_CUDA
mris_inflate_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS)
mris_inflate_cuda_LDFLAGS = $(OS_LDFLAGS) 
mris_inflate_cuda_LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(mris_inflate_cuda_CFLAGS) \
        $(CFLAGS) $(mris_inflate_cuda_LDFLAGS) $(LDFLAGS) -o $@

endif

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
