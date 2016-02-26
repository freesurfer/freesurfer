SUBDIRS =

AM_CFLAGS=-I$(top_srcdir)/include $(BOOST_CFLAGS)

AM_CXXFLAGS=-I$(top_srcdir)/include $(BOOST_CFLAGS) -Wall -Wextra

AM_LDFLAGS=

## ----
## CUDA
## ----

if HAVE_BOOST_LIBS

# BUILDCUDA is defined if configure.in finds CUDA
if BUILDCUDA

if HAVE_FERMI_GPU
NVCCFLAGS += -g -arch sm_20 --ptxas-options=-v

CUDA_CFLAGS += -DGCAMORPH_ON_GPU
endif

# rules for building cuda files
.cu.o:
	$(NVCC) -o $@ -c $< $(NVCCFLAGS) $(CUDA_CFLAGS) $(AM_CFLAGS) $(MNI_CFLAGS)


bin_PROGRAMS = test_gcamcomputelabels_boost


test_gcamcomputelabels_boost_SOURCES = test.cpp
test_gcamcomputelabels_boost_CXXFLAGS = -Wno-unused-parameter $(BOOST_CFLAGS) $(AM_CXXFLAGS)
test_gcamcomputelabels_boost_LDADD=  $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(BOOST_PO_LIB)
test_gcamcomputelabels_boost_LDFLAGS=$(OS_LDFLAGS)



endif
# endif BUILDCUDA

endif
#endif HAVE_BOOST_LIBS


# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra

