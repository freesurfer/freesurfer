SUBDIRS = gcatest \
          volumecopy \
          mrilabeltest

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


bin_PROGRAMS = gcamCompare\
               mri2nc\
               gcam_write_test gcam_write_test_cuda\
               conv1d_test conv1d_test_cuda\
               convgaussian_test convgaussian_test_cuda\
               gcam_cmp_test gcam_cmp_test_cuda\
               gcam_lle_test gcam_lle_test_cuda\
               gcam_jacobenergy_test gcam_jacobenergy_test_cuda\
               gcam_labelenergy_test gcam_labelenergy_test_cuda\
               gcam_smoothnessenergy_test gcam_smoothnessenergy_test_cuda\
               gcam_clear_test gcam_clear_test_cuda\
               gcam_applygradient_test gcam_applygradient_test_cuda\
               gcam_undogradient_test gcam_undogradient_test_cuda\
               meanfilter_test meanfilter_test_cuda\
               gcam_smoothnessterm_test gcam_smoothnessterm_test_cuda\
               gcam_jacobterm_test gcam_jacobterm_test_cuda\
               gcam_llterm_test gcam_llterm_test_cuda\
               gcam_addstatus_test gcam_addstatus_test_cuda\
               gcam_removestatus_test gcam_removestatus_test_cuda\
               gcam_labelterm_final_test gcam_labelterm_final_test_cuda\
               gcam_labelterm_postantconsist_test gcam_labelterm_postantconsist_test_cuda\
               gcam_labelterm_copydeltas_test gcam_labelterm_copydeltas_test_cuda\
               gcam_labelterm_removeoutliers_test gcam_labelterm_removeoutliers_test_cuda\
               gcam_labelterm_mainloop_test gcam_labelterm_mainloop_test_cuda\
               gcam_labelterm_test gcam_labelterm_test_cuda\
               gcam_smoothgradient_test gcam_smoothgradient_test_cuda\
               gcam_removesingularities_test gcam_removesingularities_test_cuda


gcamCompare_SOURCES = gcamCompare.cpp kahan.hpp
gcamCompare_LDADD=  $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(BOOST_PO_LIB)
gcamCompare_LDFLAGS=$(OS_LDFLAGS)

# --

gcam_write_test_SOURCES = gcam_write_test.cpp
gcam_write_test_LDADD=  $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(BOOST_PO_LIB)
gcam_write_test_LDFLAGS=$(OS_LDFLAGS)

gcam_write_test_cuda_SOURCES = gcam_write_test.cpp
gcam_write_test_cuda_CFLAGS = $(CUDA_CFLAGS) $(AM_CFLAGS) -DFS_CUDA
gcam_write_test_cuda_CXXFLAGS = $(CUDA_CXXFLAGS) $(AM_CFLAGS) -DFS_CUDA
gcam_write_test_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS) $(BOOST_PO_LIB)
gcam_write_test_cuda_LDFLAGS = $(OS_LDFLAGS) 
gcam_write_test_cuda_LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(gcam_write_test_cuda_CFLAGS) \
        $(CFLAGS) $(gcam_write_test_cuda_LDFLAGS) $(LDFLAGS) -o $@


# --

mri2nc_SOURCES = mri2nc.cpp
mri2nc_LDADD=  $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(BOOST_PO_LIB)
mri2nc_LDFLAGS=$(OS_LDFLAGS)

# --

conv1d_test_SOURCES = conv1d_test.cpp
conv1d_test_LDADD=  $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(BOOST_PO_LIB)
conv1d_test_LDFLAGS=$(OS_LDFLAGS)


conv1d_test_cuda_SOURCES = conv1d_test.cpp
conv1d_test_cuda_CFLAGS = $(CUDA_CFLAGS) $(AM_CFLAGS) -DFS_CUDA
conv1d_test_cuda_CXXFLAGS = $(CUDA_CXXFLAGS) $(AM_CFLAGS) -DFS_CUDA
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
convgaussian_test_cuda_CFLAGS = $(CUDA_CFLAGS) $(AM_CFLAGS) -DFS_CUDA
convgaussian_test_cuda_CXXFLAGS = $(CUDA_CXXFLAGS) $(AM_CFLAGS) -DFS_CUDA
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
gcam_cmp_test_cuda_CFLAGS = $(CUDA_CFLAGS)  $(AM_CFLAGS) -DFS_CUDA
gcam_cmp_test_cuda_CXXFLAGS = $(CUDA_CXXFLAGS) $(AM_CFLAGS) -DFS_CUDA
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
gcam_lle_test_cuda_CFLAGS = $(CUDA_CFLAGS)  $(AM_CFLAGS) -DFS_CUDA
gcam_lle_test_cuda_CXXFLAGS = $(CUDA_CXXFLAGS) $(AM_CFLAGS) -DFS_CUDA
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
gcam_jacobenergy_test_cuda_CFLAGS = $(CUDA_CFLAGS)  $(AM_CFLAGS) -DFS_CUDA
gcam_jacobenergy_test_cuda_CXXFLAGS = $(CUDA_CXXFLAGS) $(AM_CFLAGS) -DFS_CUDA
gcam_jacobenergy_test_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS) $(BOOST_PO_LIB)
gcam_jacobenergy_test_cuda_LDFLAGS = $(OS_LDFLAGS) 
gcam_jacobenergy_test_cuda_LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(gcam_jacobenergy_test_cuda_CFLAGS) \
        $(CFLAGS) $(gcam_jacobenergy_test_cuda_LDFLAGS) $(LDFLAGS) -o $@

# --

gcam_labelenergy_test_SOURCES = gcam_labelenergy_test.cpp
gcam_labelenergy_test_LDADD=  $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(BOOST_PO_LIB)
gcam_labelenergy_test_LDFLAGS=$(OS_LDFLAGS)


gcam_labelenergy_test_cuda_SOURCES = gcam_labelenergy_test.cpp
gcam_labelenergy_test_cuda_CFLAGS = $(CUDA_CFLAGS) $(AM_CFLAGS) -DFS_CUDA
gcam_labelenergy_test_cuda_CXXFLAGS = $(CUDA_CXXFLAGS) $(AM_CFLAGS) -DFS_CUDA
gcam_labelenergy_test_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS) $(BOOST_PO_LIB)
gcam_labelenergy_test_cuda_LDFLAGS = $(OS_LDFLAGS) 
gcam_labelenergy_test_cuda_LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(gcam_labelenergy_test_cuda_CFLAGS) \
        $(CFLAGS) $(gcam_labelenergy_test_cuda_LDFLAGS) $(LDFLAGS) -o $@


# --

gcam_smoothnessenergy_test_SOURCES = gcam_smoothnessenergy_test.cpp
gcam_smoothnessenergy_test_LDADD=  $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(BOOST_PO_LIB)
gcam_smoothnessenergy_test_LDFLAGS=$(OS_LDFLAGS)


gcam_smoothnessenergy_test_cuda_SOURCES = gcam_smoothnessenergy_test.cpp
gcam_smoothnessenergy_test_cuda_CFLAGS = $(CUDA_CFLAGS)  $(AM_CFLAGS) -DFS_CUDA
gcam_smoothnessenergy_test_cuda_CXXFLAGS = $(CUDA_CXXFLAGS) $(AM_CFLAGS) -DFS_CUDA
gcam_smoothnessenergy_test_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS) $(BOOST_PO_LIB)
gcam_smoothnessenergy_test_cuda_LDFLAGS = $(OS_LDFLAGS) 
gcam_smoothnessenergy_test_cuda_LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(gcam_smoothnessenergy_test_cuda_CFLAGS) \
        $(CFLAGS) $(gcam_smoothnessenergy_test_cuda_LDFLAGS) $(LDFLAGS) -o $@


# --

gcam_clear_test_SOURCES = gcam_clear_test.cpp
gcam_clear_test_LDADD=  $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(BOOST_PO_LIB)
gcam_clear_test_LDFLAGS=$(OS_LDFLAGS)


gcam_clear_test_cuda_SOURCES = gcam_clear_test.cpp
gcam_clear_test_cuda_CFLAGS = $(CUDA_CFLAGS) $(AM_CFLAGS) -DFS_CUDA
gcam_clear_test_cuda_CXXFLAGS = $(CUDA_CXXFLAGS) $(AM_CFLAGS) -DFS_CUDA
gcam_clear_test_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS) $(BOOST_PO_LIB)
gcam_clear_test_cuda_LDFLAGS = $(OS_LDFLAGS) 
gcam_clear_test_cuda_LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(gcam_clear_test_cuda_CFLAGS) \
        $(CFLAGS) $(gcam_clear_test_cuda_LDFLAGS) $(LDFLAGS) -o $@


# --

gcam_applygradient_test_SOURCES = gcam_applygradient_test.cpp
gcam_applygradient_test_LDADD=  $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(BOOST_PO_LIB)
gcam_applygradient_test_LDFLAGS=$(OS_LDFLAGS)


gcam_applygradient_test_cuda_SOURCES = gcam_applygradient_test.cpp
gcam_applygradient_test_cuda_CFLAGS = $(CUDA_CFLAGS) $(AM_CFLAGS) -DFS_CUDA
gcam_applygradient_test_cuda_CXXFLAGS = $(CUDA_CXXFLAGS) $(AM_CFLAGS) -DFS_CUDA
gcam_applygradient_test_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS) $(BOOST_PO_LIB)
gcam_applygradient_test_cuda_LDFLAGS = $(OS_LDFLAGS) 
gcam_applygradient_test_cuda_LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(gcam_applygradient_test_cuda_CFLAGS) \
        $(CFLAGS) $(gcam_applygradient_test_cuda_LDFLAGS) $(LDFLAGS) -o $@

# --

gcam_undogradient_test_SOURCES = gcam_undogradient_test.cpp
gcam_undogradient_test_LDADD=  $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(BOOST_PO_LIB)
gcam_undogradient_test_LDFLAGS=$(OS_LDFLAGS)


gcam_undogradient_test_cuda_SOURCES = gcam_undogradient_test.cpp
gcam_undogradient_test_cuda_CFLAGS = $(CUDA_CFLAGS) $(AM_CFLAGS) -DFS_CUDA
gcam_undogradient_test_cuda_CXXFLAGS = $(CUDA_CXXFLAGS) $(AM_CFLAGS) -DFS_CUDA
gcam_undogradient_test_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS) $(BOOST_PO_LIB)
gcam_undogradient_test_cuda_LDFLAGS = $(OS_LDFLAGS) 
gcam_undogradient_test_cuda_LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(gcam_undogradient_test_cuda_CFLAGS) \
        $(CFLAGS) $(gcam_undogradient_test_cuda_LDFLAGS) $(LDFLAGS) -o $@

# --

meanfilter_test_SOURCES=meanfilter_test.cpp
meanfilter_test_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(BOOST_PO_LIB)
meanfilter_test_LDFLAGS=$(OS_LDFLAGS)

meanfilter_test_cuda_SOURCES = meanfilter_test.cpp
meanfilter_test_cuda_CFLAGS = $(AM_CFLAGS) $(CUDA_CFLAGS) -DFS_CUDA
meanfilter_test_cuda_CXXFLAGS = $(CUDA_CXXFLAGS) $(AM_CFLAGS) -DFS_CUDA
meanfilter_test_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS) $(BOOST_PO_LIB)
meanfilter_test_cuda_LDFLAGS = $(OS_LDFLAGS) 
meanfilter_test_cuda_LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(meanfilter_test_cuda_CFLAGS) \
        $(CFLAGS) $(meanfilter_test_cuda_LDFLAGS) $(LDFLAGS) -o $@


# --

gcam_smoothnessterm_test_SOURCES = gcam_smoothnessterm_test.cpp
gcam_smoothnessterm_test_LDADD=  $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(BOOST_PO_LIB)
gcam_smoothnessterm_test_LDFLAGS=$(OS_LDFLAGS)


gcam_smoothnessterm_test_cuda_SOURCES = gcam_smoothnessterm_test.cpp
gcam_smoothnessterm_test_cuda_CFLAGS = $(CUDA_CFLAGS) $(AM_CFLAGS) -DFS_CUDA
gcam_smoothnessterm_test_cuda_CXXFLAGS = $(CUDA_CXXFLAGS) $(AM_CFLAGS) -DFS_CUDA
gcam_smoothnessterm_test_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS) $(BOOST_PO_LIB)
gcam_smoothnessterm_test_cuda_LDFLAGS = $(OS_LDFLAGS) 
gcam_smoothnessterm_test_cuda_LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(gcam_smoothnessterm_test_cuda_CFLAGS) \
        $(CFLAGS) $(gcam_smoothnessterm_test_cuda_LDFLAGS) $(LDFLAGS) -o $@

# --

gcam_jacobterm_test_SOURCES = gcam_jacobterm_test.cpp
gcam_jacobterm_test_LDADD=  $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(BOOST_PO_LIB)
gcam_jacobterm_test_LDFLAGS=$(OS_LDFLAGS)


gcam_jacobterm_test_cuda_SOURCES = gcam_jacobterm_test.cpp
gcam_jacobterm_test_cuda_CFLAGS = $(CUDA_CFLAGS) $(AM_CFLAGS) -DFS_CUDA
gcam_jacobterm_test_cuda_CXXFLAGS = $(CUDA_CXXFLAGS) $(AM_CFLAGS) -DFS_CUDA
gcam_jacobterm_test_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS) $(BOOST_PO_LIB)
gcam_jacobterm_test_cuda_LDFLAGS = $(OS_LDFLAGS) 
gcam_jacobterm_test_cuda_LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(gcam_jacobterm_test_cuda_CFLAGS) \
        $(CFLAGS) $(gcam_jacobterm_test_cuda_LDFLAGS) $(LDFLAGS) -o $@

# --

gcam_llterm_test_SOURCES = gcam_llterm_test.cpp
gcam_llterm_test_LDADD=  $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(BOOST_PO_LIB)
gcam_llterm_test_LDFLAGS=$(OS_LDFLAGS)


gcam_llterm_test_cuda_SOURCES = gcam_llterm_test.cpp
gcam_llterm_test_cuda_CFLAGS = $(CUDA_CFLAGS) $(AM_CFLAGS) -DFS_CUDA
gcam_llterm_test_cuda_CXXFLAGS = $(CUDA_CXXFLAGS) $(AM_CFLAGS) -DFS_CUDA
gcam_llterm_test_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS) $(BOOST_PO_LIB)
gcam_llterm_test_cuda_LDFLAGS = $(OS_LDFLAGS) 
gcam_llterm_test_cuda_LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(gcam_llterm_test_cuda_CFLAGS) \
        $(CFLAGS) $(gcam_llterm_test_cuda_LDFLAGS) $(LDFLAGS) -o $@

# --

gcam_labelterm_final_test_SOURCES = gcam_labelterm_final_test.cpp
gcam_labelterm_final_test_LDADD=  $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(BOOST_PO_LIB)
gcam_labelterm_final_test_LDFLAGS=$(OS_LDFLAGS)


gcam_labelterm_final_test_cuda_SOURCES = gcam_labelterm_final_test.cpp
gcam_labelterm_final_test_cuda_CFLAGS = $(CUDA_CFLAGS) $(AM_CFLAGS) -DFS_CUDA
gcam_labelterm_final_test_cuda_CXXFLAGS = $(CUDA_CXXFLAGS) $(AM_CFLAGS) -DFS_CUDA
gcam_labelterm_final_test_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS) $(BOOST_PO_LIB)
gcam_labelterm_final_test_cuda_LDFLAGS = $(OS_LDFLAGS) 
gcam_labelterm_final_test_cuda_LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(gcam_labelterm_final_test_cuda_CFLAGS) \
        $(CFLAGS) $(gcam_labelterm_final_test_cuda_LDFLAGS) $(LDFLAGS) -o $@


# --

gcam_labelterm_postantconsist_test_SOURCES = gcam_labelterm_postantconsist_test.cpp
gcam_labelterm_postantconsist_test_LDADD=  $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(BOOST_PO_LIB)
gcam_labelterm_postantconsist_test_LDFLAGS=$(OS_LDFLAGS)


gcam_labelterm_postantconsist_test_cuda_SOURCES = gcam_labelterm_postantconsist_test.cpp
gcam_labelterm_postantconsist_test_cuda_CFLAGS = $(CUDA_CFLAGS) $(AM_CFLAGS) -DFS_CUDA
gcam_labelterm_postantconsist_test_cuda_CXXFLAGS = $(CUDA_CXXFLAGS) $(AM_CFLAGS) -DFS_CUDA
gcam_labelterm_postantconsist_test_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS) $(BOOST_PO_LIB)
gcam_labelterm_postantconsist_test_cuda_LDFLAGS = $(OS_LDFLAGS) 
gcam_labelterm_postantconsist_test_cuda_LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(gcam_labelterm_postantconsist_test_cuda_CFLAGS) \
        $(CFLAGS) $(gcam_labelterm_postantconsist_test_cuda_LDFLAGS) $(LDFLAGS) -o $@

# --

gcam_labelterm_copydeltas_test_SOURCES = gcam_labelterm_copydeltas_test.cpp
gcam_labelterm_copydeltas_test_LDADD=  $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(BOOST_PO_LIB)
gcam_labelterm_copydeltas_test_LDFLAGS=$(OS_LDFLAGS)


gcam_labelterm_copydeltas_test_cuda_SOURCES = gcam_labelterm_copydeltas_test.cpp
gcam_labelterm_copydeltas_test_cuda_CFLAGS = $(CUDA_CFLAGS) $(AM_CFLAGS) -DFS_CUDA
gcam_labelterm_copydeltas_test_cuda_CXXFLAGS = $(CUDA_CXXFLAGS) $(AM_CFLAGS) -DFS_CUDA
gcam_labelterm_copydeltas_test_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS) $(BOOST_PO_LIB)
gcam_labelterm_copydeltas_test_cuda_LDFLAGS = $(OS_LDFLAGS) 
gcam_labelterm_copydeltas_test_cuda_LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(gcam_labelterm_copydeltas_test_cuda_CFLAGS) \
        $(CFLAGS) $(gcam_labelterm_copydeltas_test_cuda_LDFLAGS) $(LDFLAGS) -o $@


# --

gcam_labelterm_removeoutliers_test_SOURCES = gcam_labelterm_removeoutliers_test.cpp
gcam_labelterm_removeoutliers_test_LDADD=  $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(BOOST_PO_LIB)
gcam_labelterm_removeoutliers_test_LDFLAGS=$(OS_LDFLAGS)


gcam_labelterm_removeoutliers_test_cuda_SOURCES = gcam_labelterm_removeoutliers_test.cpp
gcam_labelterm_removeoutliers_test_cuda_CFLAGS = $(CUDA_CFLAGS) $(AM_CFLAGS) -DFS_CUDA
gcam_labelterm_removeoutliers_test_cuda_CXXFLAGS = $(CUDA_CXXFLAGS) $(AM_CFLAGS) -DFS_CUDA
gcam_labelterm_removeoutliers_test_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS) $(BOOST_PO_LIB)
gcam_labelterm_removeoutliers_test_cuda_LDFLAGS = $(OS_LDFLAGS) 
gcam_labelterm_removeoutliers_test_cuda_LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(gcam_labelterm_removeoutliers_test_cuda_CFLAGS) \
        $(CFLAGS) $(gcam_labelterm_removeoutliers_test_cuda_LDFLAGS) $(LDFLAGS) -o $@


# --

gcam_labelterm_mainloop_test_SOURCES = gcam_labelterm_mainloop_test.cpp
gcam_labelterm_mainloop_test_LDADD=  $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(BOOST_PO_LIB)
gcam_labelterm_mainloop_test_LDFLAGS=$(OS_LDFLAGS)


gcam_labelterm_mainloop_test_cuda_SOURCES = gcam_labelterm_mainloop_test.cpp
gcam_labelterm_mainloop_test_cuda_CFLAGS = $(CUDA_CFLAGS) $(AM_CFLAGS) -DFS_CUDA
gcam_labelterm_mainloop_test_cuda_CXXFLAGS = $(CUDA_CXXFLAGS) $(AM_CFLAGS) -DFS_CUDA
gcam_labelterm_mainloop_test_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS) $(BOOST_PO_LIB)
gcam_labelterm_mainloop_test_cuda_LDFLAGS = $(OS_LDFLAGS) 
gcam_labelterm_mainloop_test_cuda_LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(gcam_labelterm_mainloop_test_cuda_CFLAGS) \
        $(CFLAGS) $(gcam_labelterm_mainloop_test_cuda_LDFLAGS) $(LDFLAGS) -o $@



# --

gcam_labelterm_test_SOURCES = gcam_labelterm_test.cpp
gcam_labelterm_test_LDADD=  $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(BOOST_PO_LIB)
gcam_labelterm_test_LDFLAGS=$(OS_LDFLAGS)


gcam_labelterm_test_cuda_SOURCES = gcam_labelterm_test.cpp
gcam_labelterm_test_cuda_CFLAGS = $(CUDA_CFLAGS) $(AM_CFLAGS) -DFS_CUDA
gcam_labelterm_test_cuda_CXXFLAGS = $(CUDA_CXXFLAGS) $(AM_CFLAGS) -DFS_CUDA
gcam_labelterm_test_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS) $(BOOST_PO_LIB)
gcam_labelterm_test_cuda_LDFLAGS = $(OS_LDFLAGS) 
gcam_labelterm_test_cuda_LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(gcam_labelterm_test_cuda_CFLAGS) \
        $(CFLAGS) $(gcam_labelterm_test_cuda_LDFLAGS) $(LDFLAGS) -o $@



# --

gcam_addstatus_test_SOURCES = gcam_addstatus_test.cpp
gcam_addstatus_test_LDADD=  $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(BOOST_PO_LIB)
gcam_addstatus_test_LDFLAGS=$(OS_LDFLAGS)


gcam_addstatus_test_cuda_SOURCES = gcam_addstatus_test.cpp
gcam_addstatus_test_cuda_CFLAGS = $(CUDA_CFLAGS) $(AM_CFLAGS) -DFS_CUDA
gcam_addstatus_test_cuda_CXXFLAGS = $(CUDA_CXXFLAGS) $(AM_CFLAGS) -DFS_CUDA
gcam_addstatus_test_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS) $(BOOST_PO_LIB)
gcam_addstatus_test_cuda_LDFLAGS = $(OS_LDFLAGS) 
gcam_addstatus_test_cuda_LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(gcam_addstatus_test_cuda_CFLAGS) \
        $(CFLAGS) $(gcam_addstatus_test_cuda_LDFLAGS) $(LDFLAGS) -o $@


# --

gcam_removestatus_test_SOURCES = gcam_removestatus_test.cpp
gcam_removestatus_test_LDADD=  $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(BOOST_PO_LIB)
gcam_removestatus_test_LDFLAGS=$(OS_LDFLAGS)


gcam_removestatus_test_cuda_SOURCES = gcam_removestatus_test.cpp
gcam_removestatus_test_cuda_CFLAGS = $(CUDA_CFLAGS) $(AM_CFLAGS) -DFS_CUDA
gcam_removestatus_test_cuda_CXXFLAGS = $(CUDA_CXXFLAGS) $(AM_CFLAGS) -DFS_CUDA
gcam_removestatus_test_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS) $(BOOST_PO_LIB)
gcam_removestatus_test_cuda_LDFLAGS = $(OS_LDFLAGS) 
gcam_removestatus_test_cuda_LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(gcam_removestatus_test_cuda_CFLAGS) \
        $(CFLAGS) $(gcam_removestatus_test_cuda_LDFLAGS) $(LDFLAGS) -o $@

# --

gcam_smoothgradient_test_SOURCES = gcam_smoothgradient_test.cpp
gcam_smoothgradient_test_LDADD=  $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(BOOST_PO_LIB)
gcam_smoothgradient_test_LDFLAGS=$(OS_LDFLAGS)


gcam_smoothgradient_test_cuda_SOURCES = gcam_smoothgradient_test.cpp
gcam_smoothgradient_test_cuda_CFLAGS = $(CUDA_CFLAGS) $(AM_CFLAGS) -DFS_CUDA
gcam_smoothgradient_test_cuda_CXXFLAGS = $(CUDA_CXXFLAGS) $(AM_CFLAGS) -DFS_CUDA
gcam_smoothgradient_test_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS) $(BOOST_PO_LIB)
gcam_smoothgradient_test_cuda_LDFLAGS = $(OS_LDFLAGS) 
gcam_smoothgradient_test_cuda_LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(gcam_smoothgradient_test_cuda_CFLAGS) \
        $(CFLAGS) $(gcam_smoothgradient_test_cuda_LDFLAGS) $(LDFLAGS) -o $@

# --

gcam_removesingularities_test_SOURCES = gcam_removesingularities_test.cpp
gcam_removesingularities_test_LDADD=  $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(BOOST_PO_LIB)
gcam_removesingularities_test_LDFLAGS=$(OS_LDFLAGS)


gcam_removesingularities_test_cuda_SOURCES = gcam_removesingularities_test.cpp
gcam_removesingularities_test_cuda_CFLAGS = $(CUDA_CFLAGS) $(AM_CFLAGS) -DFS_CUDA
gcam_removesingularities_test_cuda_CXXFLAGS = $(CUDA_CXXFLAGS) $(AM_CFLAGS) -DFS_CUDA
gcam_removesingularities_test_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS) $(BOOST_PO_LIB)
gcam_removesingularities_test_cuda_LDFLAGS = $(OS_LDFLAGS) 
gcam_removesingularities_test_cuda_LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(gcam_removesingularities_test_cuda_CFLAGS) \
        $(CFLAGS) $(gcam_removesingularities_test_cuda_LDFLAGS) $(LDFLAGS) -o $@



endif
# endif BUILDCUDA

endif
#endif HAVE_BOOST_LIBS


# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra

