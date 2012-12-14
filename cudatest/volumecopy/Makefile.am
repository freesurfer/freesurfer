AM_CFLAGS=-I$(top_srcdir)/include

AM_CXXFLAGS=-I$(top_srcdir)/include $(BOOST_CFLAGS) -Wall -Wextra

AM_LDFLAGS=

## ----
## CUDA
## ----

if HAVE_BOOST_LIBS

# BUILDCUDA is defined if configure.in finds CUDA
if BUILDCUDA

CUDA_CFLAGS += -DGCAMORPH_ON_GPU

# rules for building cuda files
.cu.o:
	$(NVCC) -o $@ -c $< $(NVCCFLAGS) $(CUDA_CFLAGS) $(AM_CFLAGS) $(MNI_CFLAGS)

bin_PROGRAMS = volume_copy_test_cuda

volume_copy_test_cuda_SOURCES = volume_copy_test.cpp\
                                volcopy.cu volcopy.hpp
volume_copy_test_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS) $(BOOST_PO_LIB)
volume_copy_test_cuda_LDFLAGS = $(OS_LDFLAGS)

endif
# endif BUILDCUDA

endif
# endif HAVE_BOOST_LIBS

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra

