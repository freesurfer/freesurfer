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


bin_PROGRAMS = gcamCompare\
               conv1d_test conv1d_test_cuda\
               convgaussian_test convgaussian_test_cuda\
               gcam_cmp_test gcam_cmp_test_cuda\
               gcam_lle_test gcam_lle_test_cuda\
               gcam_jacobenergy_test gcam_jacobenergy_test_cuda


gcamCompare_SOURCES = gcamCompare.cpp
gcamCompare_LDADD=  $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(BOOST_PO_LIB)
gcamCompare_LDFLAGS=$(OS_LDFLAGS)

# --

conv1d_test_SOURCES = conv1d_test.cpp
conv1d_test_LDADD=  $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(BOOST_PO_LIB)
conv1d_test_LDFLAGS=$(OS_LDFLAGS)


conv1d_test_cuda_SOURCES = conv1d_test.cpp
conv1d_test_cuda_CFLAGS = $(AM_CFLAGS)  $(CUDA_CFLAGS) -DFS_CUDA
conv1d_test_cuda_CXXFLAGS = $(AM_CFLAGS) $(CUDA_CFLAGS) -DFS_CUDA
conv1d_test_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS) $(BOOST_PO_LIB)
conv1d_test_cuda_LDFLAGS = $(OS_LDFLAGS) 
conv1d_test_cuda_LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(conv1d_test_cuda_CFLAGS) \
        $(CFLAGS) $(conv1d_test_cuda_LDFLAGS) $(LDFLAGS) -o $@

# --

convgaussian_test_SOURCES = convgaussian_test.cpp
convgaussian_test_LDADD=  $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(BOOST_PO_LIB)
convgaussian_test_LDFLAGS=$(OS_LDFLAGS)


convgaussian_test_cuda_SOURCES = convgaussian_test.cpp
convgaussian_test_cuda_CFLAGS = $(AM_CFLAGS)  $(CUDA_CFLAGS) -DFS_CUDA
convgaussian_test_cuda_CXXFLAGS = $(AM_CFLAGS) $(CUDA_CFLAGS) -DFS_CUDA
convgaussian_test_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS) $(BOOST_PO_LIB)
convgaussian_test_cuda_LDFLAGS = $(OS_LDFLAGS) 
convgaussian_test_cuda_LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(convgaussian_test_cuda_CFLAGS) \
        $(CFLAGS) $(convgaussian_test_cuda_LDFLAGS) $(LDFLAGS) -o $@

# --

gcam_cmp_test_SOURCES = gcam_cmp_test.cpp
gcam_cmp_test_LDADD=  $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(BOOST_PO_LIB)
gcam_cmp_test_LDFLAGS=$(OS_LDFLAGS)


gcam_cmp_test_cuda_SOURCES = gcam_cmp_test.cpp
gcam_cmp_test_cuda_CFLAGS = $(AM_CFLAGS)  $(CUDA_CFLAGS) -DFS_CUDA
gcam_cmp_test_cuda_CXXFLAGS = $(AM_CFLAGS) $(CUDA_CFLAGS) -DFS_CUDA
gcam_cmp_test_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS) $(BOOST_PO_LIB)
gcam_cmp_test_cuda_LDFLAGS = $(OS_LDFLAGS) 
gcam_cmp_test_cuda_LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(gcam_cmp_test_cuda_CFLAGS) \
        $(CFLAGS) $(gcam_cmp_test_cuda_LDFLAGS) $(LDFLAGS) -o $@

# --

gcam_lle_test_SOURCES = gcam_lle_test.cpp
gcam_lle_test_LDADD=  $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(BOOST_PO_LIB)
gcam_lle_test_LDFLAGS=$(OS_LDFLAGS)


gcam_lle_test_cuda_SOURCES = gcam_lle_test.cpp
gcam_lle_test_cuda_CFLAGS = $(AM_CFLAGS)  $(CUDA_CFLAGS) -DFS_CUDA
gcam_lle_test_cuda_CXXFLAGS = $(AM_CFLAGS) $(CUDA_CFLAGS) -DFS_CUDA
gcam_lle_test_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS) $(BOOST_PO_LIB)
gcam_lle_test_cuda_LDFLAGS = $(OS_LDFLAGS) 
gcam_lle_test_cuda_LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(gcam_lle_test_cuda_CFLAGS) \
        $(CFLAGS) $(gcam_lle_test_cuda_LDFLAGS) $(LDFLAGS) -o $@

# --

gcam_jacobenergy_test_SOURCES = gcam_jacobenergy_test.cpp
gcam_jacobenergy_test_LDADD=  $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(BOOST_PO_LIB)
gcam_jacobenergy_test_LDFLAGS=$(OS_LDFLAGS)


gcam_jacobenergy_test_cuda_SOURCES = gcam_jacobenergy_test.cpp
gcam_jacobenergy_test_cuda_CFLAGS = $(AM_CFLAGS)  $(CUDA_CFLAGS) -DFS_CUDA
gcam_jacobenergy_test_cuda_CXXFLAGS = $(AM_CFLAGS) $(CUDA_CFLAGS) -DFS_CUDA
gcam_jacobenergy_test_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS) $(BOOST_PO_LIB)
gcam_jacobenergy_test_cuda_LDFLAGS = $(OS_LDFLAGS) 
gcam_jacobenergy_test_cuda_LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(gcam_jacobenergy_test_cuda_CFLAGS) \
        $(CFLAGS) $(gcam_jacobenergy_test_cuda_LDFLAGS) $(LDFLAGS) -o $@




endif
# endif BUILDCUDA


# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra

