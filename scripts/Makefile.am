#
# Makefile.am for scripts
#
dist_bin_SCRIPTS = \
archive_pull \
aparc2feat \
aseg2feat \
aparcstats2table \
asegstats2table \
bblabel \
browse-minc-header.tcl \
browse-sessions \
check_subject \
clear_fs_env.csh \
cor_to_minc \
csurf \
dcmdir-info-mgh \
deface_subject \
epidewarp.fsl\
feat2surf \
findsession \
fixup_mni_paths \
fix_subject \
fix_subject_corrected \
fix_subject_corrected-lh \
fix_subject_corrected-rh \
fix_subject-lh \
fix_subject_on_seychelles \
fix_subject-rh \
freesurfer.csh \
fs_lib_check \
fsl_rigid_register \
fslregister \
getfullpath \
grad_unwarp \
inflate_subject \
inflate_subject3 \
inflate_subject-lh \
inflate_subject_new \
inflate_subject_new-lh \
inflate_subject_new-rh \
inflate_subject-rh \
inflate_subject_sc \
isanalyze \
isnifti \
label_child \
label_elderly_subject \
label_subject \
label_subject_flash \
label_subject_mixed \
make_average_surface \
make_average_volume \
make_average_subject \
make-segvol-table \
map_all_labels \
map_all_labels-lh \
map_central_sulcus \
minc2seqinfo \
mkfsdist \
mkima_index.tcl \
mkmnc_index.tcl \
mksubjdirs \
mksurfatlas \
morph_only_subject \
morph_only_subject-lh \
morph_only_subject-rh \
morph_rgb-lh \
morph_rgb-rh \
morph_subject \
morph_subject-lh \
morph_subject_on_seychelles \
morph_subject-rh \
morph_tables-lh \
morph_tables-rh \
mri_nu_correct.mni\
mri-func2sph \
mri-funcvits \
mri_mergelabels \
mri_motion_correct \
mri_motion_correct2 \
mri_motion_correct.fsl \
mri-sph2surf \
mris_preproc \
mris_volsmooth \
ms_refine_subject \
printdat \
progressbar.tcl \
fsr-checkxopts \
fsr-getxopts \
recon-all \
recon-all-dev \
recon-all.2005 \
reg-feat2anat \
register_child \
register.csh \
register_elderly_subject \
register_subject \
register_subject_flash \
register_subject_mixed \
reinflate_subject \
reinflate_subject-lh \
reinflate_subject-rh \
remove_talairach \
renormalize_subject \
renormalize_subject_keep_editting \
renormalize_T1_subject \
reregister_subject_mixed \
run-qdec-glm \
segment_monkey \
segment_subject \
segment_subject_notal \
segment_subject_notal2 \
segment_subject_old_skull_strip \
segment_subject_sc \
segment_subject_talmgh \
setlabelstat\
show_tal \
sphere_subject \
sphere_subject-lh \
sphere_subject-rh \
spmmat2register \
spm_t_to_b \
spmregister \
talairach \
talairach2 \
talairach_mgh \
tal_compare \
tkmedit_wrapper \
unpack_ima1.tcl \
unpackimadir \
unpackimadir2 \
unpack_ima.tcl \
unpackmincdir \
unpack_mnc.tcl \
unpacksdcmdir \
wfilemask \
xsanatreg

tcldir = $(prefix)/lib/tcl
tcl_SCRIPTS= \
curv-flat.tcl \
curv-views.tcl \
dara_data.tcl \
dipmovie-flat.tcl \
eccen-flat.tcl \
eccen-views.tcl \
fs-flat.tcl \
fs-make.tcl \
fs-views.tcl \
mk-fieldsign.tcl \
offsetmovie.tcl \
phasemovie.tcl \
polar-flat.tcl \
polar-views.tcl \
readenv.tcl \
real-flat.tcl \
real-views.tcl \
saveflat.tcl \
saveviews.tcl \
segment_subject \
setdefpatchview.tcl \
solmovie-flat.tcl \
tkanalyse.tcl \
tkm_common.tcl \
tkm_dialog.tcl \
tkmedit.tcl \
tkm_functional.tcl \
tkm_graph.tcl \
tkm_wrappers.tcl \
tkregister.tcl \
tkregister2.tcl \
tksurfer.tcl \
tkUtils.tcl \
twocond-flat.tcl \
twocond-masked-views.tcl \
twocond-views.tcl \
wrappers.tcl \
fsgdfPlot.tcl

EXTRA_DIST = csurf.in $(tcl_SCRIPTS)

# notes: 
# recon-all.2005 is the renamed recon-all-nmr, which is no longer used.
# recon-all is the official recon-all.
# recon-all-dev is the 'development' version of recon-all, which,
#               when changes are tested, get merged into recon-all.

# use , as the separator so that the directory path / can be passed
SEDFLAGS = -e 's,xEXEC_PREFIXx,$(exec_prefix),g' \
	   -e 's,xPREFIXx,$(prefix),g' \
	   -e 's,xTIXWISHx,$(TIXWISH),g'

csurf: csurf.in
	rm -f $@
	sed < $< $(SEDFLAGS) > $@
	chmod +x $@

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
