#
# Makefile.am
#
AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

check_PROGRAMS=mnitest i2rtest icotest surftest extest mghxform inftest difftool mymri_mask mris_average_on_one_surface mris_diff_on_surface mri_ms_fcm mri_compute_stats mri_ms_LDA mri_ms_LDA_local mri_fslmat_to_lta mri_concatenate_lta checkWfile mri_transform_to_COR mri_ms_LDA_local_all mri_ms_LDA_all mri_seg_concensus mri_seg_volume mri_ms_LDA_whole mri_ms_LDA_test mris_fastmarching mri_compute_all_CNR mri_ms_LDA_all_test mri_ms_LDA_one mri_ms_compute_CNR mris_thickness_diff mris_diff mri_seg_concensus2 mris_compute_thickness mri_modify_seg mri_ms_compute_lesion_CNR mri_apply_INU mri_apply_dura_mask mri_ms_EM mri_edit_wm mri_seg_concensus_tmp mri_apply_EM_mask mri_apply_aseg_mask mri_apply_aseg_mask2 mri_apply_EM_mask2 mri_extract

TESTS=mnitest i2rtest icotest extest mghxform inftest checkanalyze

mymri_mask_SOURCE=mymri_mask.c
mris_average_on_one_surface_SOURCE=mris_average_on_one_surface.c
mris_diff_on_surface_SOURCE=mris_diff_on_surface.c
# mris_average_on_one_surface_LDFLAGS=-all-static
mri_ms_fcm_SOURCE=mri_ms_fcm.c
mri_edit_wm_SOURCE=mri_edit_wm.c
mri_apply_EM_mask_SOURCE=mri_apply_EM_mask.c
mri_ms_EM_SOURCE=mri_ms_EM.c
mri_ms_LDA_SOURCE=mri_ms_LDA.c
mris_fastmarching_SOURCES=mris_fastmarching.c fmarchmesh.c heap.c
mri_ms_LDA_all_SOURCE=mri_ms_LDA_all.c
mri_ms_LDA_all_test_SOURCE=mri_ms_LDA_all_test.c
mri_ms_LDA_one_SOURCE=mri_ms_LDA_one.c
mri_ms_LDA_whole_SOURCE=mri_ms_LDA_whole.c
mri_ms_LDA_local_SOURCE=mri_ms_LDA_local.c
mri_ms_LDA_local_all_SOURCE=mri_ms_LDA_local_all.c
mri_ms_LDA_test_SOURCE=mri_ms_LDA_test.c
mri_apply_INU_SOURCE=mri_apply_INU.c
mri_apply_dura_mask_SOURCE=mri_apply_dura_mask.c
mri_ms_compute_CNR_SOURCE=mri_ms_compute_CNR.c
mri_ms_compute_lesion_CNR_SOURCE=mri_ms_compute_lesion_CNR.c
checkWfile_SOURCE=checkWfile.c
mri_compute_stats_SOURCE=mri_compute_stats.c
mri_compute_all_CNR_SOURCE=mri_compute_all_CNR.c
mri_fslmat_to_lta_SOURCE=mri_fslmat_to_lta.c
mri_transform_to_COR_SOURCE=mri_transform_to_COR.c
mri_concatenate_lta_SOURCE=mri_concatenate_lta.c
mri_seg_concensus_SOURCE=mri_seg_concensus.c
mri_modify_seg_SOURCE=mri_modify_seg.c
mri_seg_concensus2_SOURCE=mri_seg_concensus2.c
mri_seg_concensus_tmp_SOURCE=mri_seg_concensus_tmp.c
mri_seg_volume_SOURCE=mri_seg_volume.c
mris_thickness_diff_SOURCES=mris_thickness_diff.c
mris_diff_SOURCES=mris_diff.cpp fmarching3dnband.c heap.c
mri_ms_compute_CNR_SOURCES=mri_ms_compute_CNR.c
mris_compute_thickness_SOURCES=mris_compute_thickness.cpp
mri_apply_aseg_mask_SOURCE = mri_apply_aseg_mask.c
mri_apply_aseg_mask2_SOURCE = mri_apply_aseg_mask2.c
mri_apply_EM_mask2_SOURCE = mri_apply_EM_mask2.c
mri_extract_SOURCE=mri_extract.c
mnitest_SOURCES=mnitest.cpp
i2rtest_SOURCES=i2rtest.cpp
icotest_SOURCES=icotest.cpp
surftest_SOURCES=surftest.cpp
extest_SOURCES=extest.cpp
mghxform_SOURCES=mghxform.cpp
inftest_SOURCES=inftest.cpp
difftool_SOURCES=difftool.cpp

LDADD=$(addprefix $(top_srcdir)/, $(LIBS_MGH))

mris_diff_LDFLAGS=-L/autofs/homes/003/xhan/lnxlib -L$(top_srcdir)/lib
mris_diff_LDADD=-lANN $(addprefix $(top_srcdir)/, $(LIBS_MGH))

mris_thickness_diff_LDFLAGS=$(OS_LDFLAGS)

mris_compute_thickness_LDFLAGS=-L/autofs/homes/003/xhan/lnxlib -L$(top_srcdir)/lib
mris_compute_thickness_LDADD=-lANN $(addprefix $(top_srcdir)/, $(LIBS_MGH))

