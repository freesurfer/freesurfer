#
# Makefile.am
#
AM_CPPFLAGS=-I$(top_srcdir)/include \
	-I$(top_srcdir)/include/dicom \
	-I$(top_srcdir)/include/NrrdIO 

check_PROGRAMS=mnitest i2rtest icotest surftest extest \
	mghxform inftest difftool mri_compute_stats \
	mri_ms_LDA mri_fslmat_to_lta mri_concatenate_lta \
	mri_transform_to_COR  mri_compute_all_CNR \
	mri_ms_compute_CNR \
	mri_apply_INU mri_edit_wm mri_apply_EM_mask \
	mri_extract mri_ms_gca_EM checkanalyze \
	test_mri_identify test_mriio test_c_nr_wrapper \
	sc_test
#test_nr_wrapper test_matrix 

TESTS=sc_test \
	mnitest i2rtest icotest extest mghxform \
	inftest checkanalyze test_mri_identify \
	test_c_nr_wrapper
#test_mriio test_matrix


mri_edit_wm_SOURCE=mri_edit_wm.c
mri_apply_EM_mask_SOURCE=mri_apply_EM_mask.c
mri_ms_gca_EM_SOURCE=mri_ms_EM.c
mri_ms_LDA_SOURCE=mri_ms_LDA.c
mri_apply_INU_SOURCE=mri_apply_INU.c
mri_compute_stats_SOURCE=mri_compute_stats.c
mri_compute_all_CNR_SOURCE=mri_compute_all_CNR.c
mri_fslmat_to_lta_SOURCE=mri_fslmat_to_lta.c
mri_transform_to_COR_SOURCE=mri_transform_to_COR.c
mri_concatenate_lta_SOURCE=mri_concatenate_lta.c
mri_ms_compute_CNR_SOURCES=mri_ms_compute_CNR.c
#mris_compute_thickness_SOURCES=mris_compute_thickness.cpp
mri_extract_SOURCE=mri_extract.c
mnitest_SOURCES=mnitest.cpp
i2rtest_SOURCES=i2rtest.cpp
icotest_SOURCES=icotest.cpp
surftest_SOURCES=surftest.cpp
extest_SOURCES=extest.cpp
mghxform_SOURCES=mghxform.cpp
inftest_SOURCES=inftest.cpp
difftool_SOURCES=difftool.cpp
checkanalyze_SOURCES=checkanalyze.cpp
test_mri_identify_SOURCES=test_mri_identify.cpp
#test_nr_wrapper_SOURCES=../nr_wrapper_open_source.cpp test_nr_wrapper.cpp
test_c_nr_wrapper_SOURCES=../nr_wrapper_open_source.cpp test_c_nr_wrapper.c
test_mriio_SOURCES=test_mriio.cpp
#test_matrix_SOURCES=../nr_wrapper_open_source.cpp test_matrix.cpp
sc_test_SOURCES=sc_test.c

#mris_thickness_diff_LDFLAGS=$(OS_LDFLAGS)

#mris_compute_thickness_LDFLAGS=-L/autofs/homes/003/xhan/lnxlib
#mris_compute_thickness_LDADD=-lANN $(addprefix $(top_srcdir)/, $(LIBS_MGH))

#test_nr_wrapper_LDFLAGS=-ldl
#test_nr_wrapper_LDADD=$(addprefix $(top_srcdir)/, $(LIBS_MGH) ) \
#	-lcppunit -lvnl_algo -lvnl -lnetlib

#test_matrix_LDFLAGS=-ldl
#test_matrix_LDADD=$(addprefix $(top_srcdir)/, $(LIBS_MGH) ) \
#	-lcppunit -lvnl_algo -lvnl -lnetlib

test_c_nr_wrapper_LDADD=$(addprefix $(top_srcdir)/, $(LIBS_MGH) ) \
  -lvnl_algo -lvnl -lnetlib

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES="sc_test.c"
include $(top_srcdir)/Makefile.extra
