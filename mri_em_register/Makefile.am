## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom 

AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom\
             $(BOOST_CFLAGS)\
             -Wall -Wextra


AM_LDFLAGS=

bin_PROGRAMS = mri_em_register
mri_em_register_SOURCES=mri_em_register.c
mri_em_register_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_em_register_LDFLAGS=$(OS_LDFLAGS)



## ----
## CUDA
## ----

# BUILDCUDA is defined if configure.in finds CUDA
if BUILDCUDA
# rules for building cuda files
.cu.o:
	$(NVCC) -o $@ -c $< $(NVCCFLAGS) $(AM_CFLAGS) $(MNI_CFLAGS)

bin_PROGRAMS += mri_em_register_cuda \
               meanfilter_test \
               meanfilter_test_cuda


mri_em_register_cuda_SOURCES = mri_em_register.c \
                               computelogsampleprob.cu\
                               findoptimaltranslation.cu
mri_em_register_cuda_CFLAGS = $(AM_CFLAGS) $(CUDA_CFLAGS) -DFS_CUDA
mri_em_register_cuda_CXXFLAGS = $(AM_CFLAGS) $(CUDA_CFLAGS) -DFS_CUDA
mri_em_register_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS)
mri_em_register_cuda_LDFLAGS = $(OS_LDFLAGS) 
mri_em_register_cuda_LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(mri_em_register_cuda_CFLAGS) \
        $(CFLAGS) $(mri_em_register_cuda_LDFLAGS) $(LDFLAGS) -o $@



meanfilter_test_SOURCES=meanfilter_test.cpp
meanfilter_test_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(BOOST_PO_LIB)
meanfilter_test_LDFLAGS=$(OS_LDFLAGS)

meanfilter_test_cuda_SOURCES = meanfilter_test.cpp
meanfilter_test_cuda_CFLAGS = $(AM_CFLAGS) $(CUDA_CFLAGS) -DFS_CUDA
meanfilter_test_cuda_CXXFLAGS = $(AM_CFLAGS) $(CUDA_CFLAGS) -DFS_CUDA
meanfilter_test_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS) $(BOOST_PO_LIB)
meanfilter_test_cuda_LDFLAGS = $(OS_LDFLAGS) 
meanfilter_test_cuda_LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(meanfilter_test_cuda_CFLAGS) \
        $(CFLAGS) $(meanfilter_test_cuda_LDFLAGS) $(LDFLAGS) -o $@


endif

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
