#
# Makefile.am
#

AM_CPPFLAGS=-I$(top_srcdir)/include \
	-I$(top_srcdir)/include/dicom \
	-I$(top_srcdir)/include/NrrdIO
LDADD=$(addprefix $(top_srcdir)/, $(LIBS_MGH))

CHECKS=mnitest i2rtest icotest extest \
	mghxform inftest checkanalyze \
	test_mri_identify test_c_nr_wrapper \
	sc_test

BROKEN=difftool test_mriio mri_compute_stats \
  surftest mri_ms_LDA \
	mri_fslmat_to_lta mri_concatenate_lta \
	mri_transform_to_COR mri_compute_all_CNR \
	mri_ms_compute_CNR mris_compute_thickness \
	mri_apply_INU mri_edit_wm mri_apply_EM_mask \
	mri_extract mri_ms_gca_EM

if HAVE_CPPUNIT
CPPUNIT_CHECKS=test_nr_wrapper test_matrix
check_PROGRAMS=$(CHECKS) $(CPPUNIT_CHECKS)
else
CPPUNIT_CHECKS=
check_PROGRAMS=$(CHECKS)
endif

TESTS=setup_test_data $(check_PROGRAMS) cleanup_test_data

mnitest_SOURCES=mnitest.cpp
i2rtest_SOURCES=i2rtest.cpp
icotest_SOURCES=icotest.cpp
extest_SOURCES=extest.cpp
mghxform_SOURCES=mghxform.cpp
inftest_SOURCES=inftest.cpp
checkanalyze_SOURCES=checkanalyze.cpp
test_mri_identify_SOURCES=test_mri_identify.cpp
test_c_nr_wrapper_SOURCES=test_c_nr_wrapper.c
sc_test_SOURCES=sc_test.c
#test_mriio_SOURCES=test_mriio.cpp
#surftest_SOURCES=surftest.cpp
#difftool_SOURCES=difftool.cpp
#mris_compute_thickness_SOURCES=mris_compute_thickness.cpp
#mri_edit_wm_SOURCE=mri_edit_wm.c
#mri_apply_EM_mask_SOURCE=mri_apply_EM_mask.c
#mri_ms_gca_EM_SOURCE=mri_ms_EM.c
#mri_ms_LDA_SOURCE=mri_ms_LDA.c
#mri_apply_INU_SOURCE=mri_apply_INU.c
#mri_compute_stats_SOURCE=mri_compute_stats.c
#mri_compute_all_CNR_SOURCE=mri_compute_all_CNR.c
#mri_fslmat_to_lta_SOURCE=mri_fslmat_to_lta.c
#mri_transform_to_COR_SOURCE=mri_transform_to_COR.c
#mri_concatenate_lta_SOURCE=mri_concatenate_lta.c
#mri_ms_compute_CNR_SOURCES=mri_ms_compute_CNR.c
#mri_extract_SOURCE=mri_extract.c

if HAVE_CPPUNIT

# -Wno-non-virtual-dtor suppresses valid scuba and cppunit warnings:
# "has virtual functions but non-virtual destructor"
# these warnings appear with gcc 4.#

test_nr_wrapper_SOURCES=test_nr_wrapper.cpp
test_nr_wrapper_CXXFLAGS=$(CPPUNIT_CFLAGS) -Wno-non-virtual-dtor
test_nr_wrapper_LDADD=$(addprefix $(top_srcdir)/, $(LIBS_MGH)) \
	$(CPPUNIT_LIBS) -lcppunit -ldl

test_matrix_SOURCES=test_matrix.cpp
test_matrix_CXXFLAGS=$(CPPUNIT_CFLAGS) -Wno-non-virtual-dtor
test_matrix_LDADD=$(addprefix $(top_srcdir)/, $(LIBS_MGH)) \
	$(CPPUNIT_LIBS) -lcppunit -ldl
endif

EXTRA_DIST=setup_test_data test_data.tar.gz cleanup_test_data

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
