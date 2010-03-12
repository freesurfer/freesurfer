AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom\
             $(BOOST_CFLAGS)\
             -Wall -Wextra

AM_LDFLAGS=

## ----
## CUDA
## ----

# BUILDCUDA is defined if configure.in finds CUDA
if BUILDCUDA

NVCCFLAGS += -g -arch sm_13 --ptxas-options=-v

# rules for building cuda files
.cu.o:
	$(NVCC) -o $@ -c $< $(NVCCFLAGS) $(AM_CFLAGS) $(MNI_CFLAGS)


bin_PROGRAMS = gcamCompare


gcamCompare_SOURCES = gcamCompare.cpp
gcamCompare_LDADD=  $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(BOOST_PO_LIB)
gcamCompare_LDFLAGS=$(OS_LDFLAGS)

endif
# endif BUILDCUDA


# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra

