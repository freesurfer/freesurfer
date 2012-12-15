AM_CFLAGS=-I$(top_srcdir)/include

AM_CXXFLAGS=-I$(top_srcdir)/include\
             $(BOOST_CFLAGS)\
             -Wall -Wextra

AM_LDFLAGS=

## ----
## CUDA
## ----

if HAVE_BOOST_LIBS

# BUILDCUDA is defined if configure.in finds CUDA
if BUILDCUDA

NVCCFLAGS += -g -arch sm_20 --ptxas-options=-v

CUDA_CFLAGS += -DGCAMORPH_ON_GPU

# rules for building cuda files
.cu.o:
	$(NVCC) -o $@ -c $< $(NVCCFLAGS) $(CUDA_CFLAGS) $(AM_CFLAGS) $(MNI_CFLAGS)


bin_PROGRAMS = mriMarkLabelBorderVoxelTest mriMarkLabelBorderVoxelTest_cuda\
               mriVoxInLabelPartialVolumeTest mriVoxInLabelPartialVolumeTest_cuda

# ---

mriMarkLabelBorderVoxelTest_SOURCES = mriMarkLabelBorderVoxelTest.cpp
mriMarkLabelBorderVoxelTest_LDADD = $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(BOOST_PO_LIB)
mriMarkLabelBorderVoxelTest_LDFLAGS = $(OS_LDFLAGS)

mriMarkLabelBorderVoxelTest_cuda_SOURCES = mriMarkLabelBorderVoxelTest.cpp
mriMarkLabelBorderVoxelTest_cuda_CFLAGS = $(CUDA_CFLAGS) $(AM_CFLAGS) $(BOOST_CFLAGS) -DFS_CUDA
mriMarkLabelBorderVoxelTest_cuda_CXXFLAGS = $(CUDA_CXXFLAGS) $(AM_CFLAGS) $(BOOST_CFLAGS) -DFS_CUDA
mriMarkLabelBorderVoxelTest_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS) $(BOOST_PO_LIB)
mriMarkLabelBorderVoxelTest_cuda_LDFLAGS = $(OS_LDFLAGS)
mriMarkLabelBorderVoxelTest_cuda_LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(mriMarkLabelBorderVoxelTest_cuda_CFLAGS) \
        $(CFLAGS) $(mriMarkLabelBorderVoxelTest_cuda_LDFLAGS) $(LDFLAGS) -o $@
# ---

# ---

mriVoxInLabelPartialVolumeTest_SOURCES = mriVoxInLabelPartialVolumeTest.cpp
mriVoxInLabelPartialVolumeTest_LDADD = $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(BOOST_PO_LIB)
mriVoxInLabelPartialVolumeTest_LDFLAGS = $(OS_LDFLAGS)

mriVoxInLabelPartialVolumeTest_cuda_SOURCES = mriVoxInLabelPartialVolumeTest.cpp
mriVoxInLabelPartialVolumeTest_cuda_CFLAGS = $(CUDA_CFLAGS) $(AM_CFLAGS) $(BOOST_CFLAGS) -DFS_CUDA
mriVoxInLabelPartialVolumeTest_cuda_CXXFLAGS = $(CUDA_CXXFLAGS) $(AM_CFLAGS) $(BOOST_CFLAGS) -DFS_CUDA
mriVoxInLabelPartialVolumeTest_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS) $(BOOST_PO_LIB)
mriVoxInLabelPartialVolumeTest_cuda_LDFLAGS = $(OS_LDFLAGS)
mriVoxInLabelPartialVolumeTest_cuda_LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(mriVoxInLabelPartialVolumeTest_cuda_CFLAGS) \
        $(CFLAGS) $(mriVoxInLabelParialVolumeTest_cuda_LDFLAGS) $(LDFLAGS) -o $@

endif
# endif BUILDCUDA

endif
# endif HAVE_BOOST_LIBS

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra

