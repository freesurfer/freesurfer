#
# Makefile.am
#
AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

check_PROGRAMS=mnitest i2rtest icotest surftest extest mghxform inftest difftool  mri_ms_fcm mri_compute_stats mri_ms_LDA mri_fslmat_to_lta mri_concatenate_lta checkWfile mri_transform_to_COR  mri_compute_all_CNR mri_ms_compute_CNR  mris_diff mris_compute_thickness mri_apply_INU mri_ms_EM mri_edit_wm mri_apply_EM_mask mri_extract mri_ms_gca_EM mris_indirect_morph mris_average mri_ms_compute_T2star mri_ms_optimize_seq refine_aseg seg_pallidum mri_compute_seg_overlap aseg_for_one mri_count_seg_vols mri_unmask

TESTS=mnitest i2rtest icotest extest mghxform inftest checkanalyze

aseg_for_one_SOURCES=aseg_for_one.c
mri_ms_fcm_SOURCE=mri_ms_fcm.c
mri_unmask_SOURCE=mri_unmask.c
mri_edit_wm_SOURCE=mri_edit_wm.c
mri_apply_EM_mask_SOURCE=mri_apply_EM_mask.c
mri_ms_EM_SOURCE=mri_ms_EM.c
mri_ms_gca_EM_SOURCE=mri_ms_EM.c
mri_ms_LDA_SOURCE=mri_ms_LDA.c
mri_count_seg_vols_SOURCES=mri_count_seg_vols.c	
mri_apply_INU_SOURCE=mri_apply_INU.c
checkWfile_SOURCE=checkWfile.c
mri_compute_stats_SOURCE=mri_compute_stats.c
mri_compute_seg_overlap_SOURCES=mri_compute_seg_overlap.c
mri_compute_all_CNR_SOURCE=mri_compute_all_CNR.c
mri_fslmat_to_lta_SOURCE=mri_fslmat_to_lta.c
mri_transform_to_COR_SOURCE=mri_transform_to_COR.c
mri_concatenate_lta_SOURCE=mri_concatenate_lta.c
mris_diff_SOURCES=mris_diff.cpp fmarching3dnband.c heap_vol.c
mri_ms_compute_CNR_SOURCES=mri_ms_compute_CNR.c
seg_pallidum_SOURCES=seg_pallidum.c
mri_ms_compute_T2star_SOURCES=mri_ms_compute_T2star.c
mris_compute_thickness_SOURCES=mris_compute_thickness.cpp
mri_ms_optimize_seq_SOURCE=mri_ms_optimize_seq.c
mri_extract_SOURCE=mri_extract.c
mnitest_SOURCES=mnitest.cpp
i2rtest_SOURCES=i2rtest.cpp
icotest_SOURCES=icotest.cpp
surftest_SOURCES=surftest.cpp
extest_SOURCES=extest.cpp
mghxform_SOURCES=mghxform.cpp
inftest_SOURCES=inftest.cpp
difftool_SOURCES=difftool.cpp
mris_indirect_morph_SOURCES=mris_indirect_morph.cpp
mris_average_SOURCES=mris_average.c
refine_aseg_SOURCES=refine_aseg.c

LDADD=$(addprefix $(top_srcdir)/, $(LIBS_MGH))

mris_diff_LDFLAGS=-L/autofs/homes/003/xhan/lnxlib -L$(top_srcdir)/lib
mris_diff_LDADD=-lANN $(addprefix $(top_srcdir)/, $(LIBS_MGH))

mris_indirect_morph_LDFLAGS=-L/autofs/homes/003/xhan/lnxlib -L$(top_srcdir)/lib
mris_indirect_morph_LDADD=-lANN $(addprefix $(top_srcdir)/, $(LIBS_MGH))

#mris_thickness_diff_LDFLAGS=$(OS_LDFLAGS)

mris_compute_thickness_LDFLAGS=-L/autofs/homes/003/xhan/lnxlib -L$(top_srcdir)/lib
mris_compute_thickness_LDADD=-lANN $(addprefix $(top_srcdir)/, $(LIBS_MGH))

