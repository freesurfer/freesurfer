## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom\
             $(BOOST_CFLAGS)\
             -Wall -Wextra

AM_LDFLAGS=

bin_PROGRAMS = mri_ca_register\
               conv1d_test

mri_ca_register_SOURCES=mri_ca_register.c
mri_ca_register_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH)) 
mri_ca_register_LDFLAGS=$(OS_LDFLAGS)

conv1d_test_SOURCES = conv1d_test.cpp
conv1d_test_LDADD=  $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(BOOST_PO_LIB)
conv1d_test_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra

## ----
## CUDA
## ----

# BUILDCUDA is defined if configure.in finds CUDA
if BUILDCUDA

NVCCFLAGS += -g -arch sm_13 --ptxas-options=-v

# rules for building cuda files
.cu.o:
	$(NVCC) -o $@ -c $< $(NVCCFLAGS) $(AM_CFLAGS) $(MNI_CFLAGS)

bin_PROGRAMS += mri_ca_register_cuda\
                conv1d_test_cuda
mri_ca_register_cuda_SOURCES = mri_ca_register.c\
                               devicemanagement.cu
mri_ca_register_cuda_CFLAGS = $(AM_CFLAGS)  $(CUDA_CFLAGS) -DFS_CUDA
mri_ca_register_cuda_CXXFLAGS = $(AM_CFLAGS) $(CUDA_CFLAGS) -DFS_CUDA
mri_ca_register_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS)
mri_ca_register_cuda_LDFLAGS = $(OS_LDFLAGS) 
mri_ca_register_cuda_LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(mri_em_register_cuda_CFLAGS) \
        $(CFLAGS) $(mri_em_register_cuda_LDFLAGS) $(LDFLAGS) -o $@

conv1d_test_cuda_SOURCES = conv1d_test.cpp\
                           mriconvolve1d_cuda.cu\
                           mri_transfer.cu
conv1d_test_cuda_CFLAGS = $(AM_CFLAGS)  $(CUDA_CFLAGS) -DFS_CUDA
conv1d_test_cuda_CXXFLAGS = $(AM_CFLAGS) $(CUDA_CFLAGS) -DFS_CUDA
conv1d_test_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS) $(BOOST_PO_LIB)
conv1d_test_cuda_LDFLAGS = $(OS_LDFLAGS) 
conv1d_test_cuda_LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(mri_em_register_cuda_CFLAGS) \
        $(CFLAGS) $(mri_em_register_cuda_LDFLAGS) $(LDFLAGS) -o $@

endif
