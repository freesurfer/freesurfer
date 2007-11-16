#
# Makefile.am for scripts
#
dist_bin_SCRIPTS = \
archive_pull \
aparc2feat \
aseg2feat \
aparcstats2table \
asegstats2table \
aparcstatsdiff \
asegstatsdiff \
batchjobs \
bblabel \
bedpostx_seychelles \
browse-minc-header.tcl \
browse-sessions \
bugr \
build_desikan_killiany_gcs.csh \
check_subject \
clear_fs_env.csh \
cor_to_minc \
cp-dicom \
dcmdir-info-mgh \
deface_subject \
dicom-rename \
dt_recon \
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
flip4fsl \
fname2stem \
fsl_sub_seychelles \
stem2fname \
freesurfer \
fs_lib_check \
fs_time \
fsl_rigid_register \
fslregister \
fsfget \
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
kww_setup \
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
mkheadsurf \
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
parc_atlas_jackknife_test \
petsc_setup \
plot_structure_stats.tcl \
printdat \
progressbar.tcl \
fsr-checkxopts \
fsr-getxopts \
rebuild_gca_atlas.csh \
recon-all \
recon-all-dev \
recon-all.2005 \
recon-all.makefile \
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
talsegprob \
tcl_setup \
test_recon-all.csh \
unpack_ima1.tcl \
unpackimadir \
unpackimadir2 \
unpack_ima.tcl \
unpackmincdir \
unpack_mnc.tcl \
unpacksdcmdir \
vno_match_check \
vtk_setup \
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

install-exec-hook:
	cp $(srcdir)/unpacksdcmdir $(DESTDIR)$(bindir)/unpacksdcmdir.tcl
	echo "#!/bin/tcsh -ef" > $(DESTDIR)$(bindir)/unpacksdcmdir
	echo "source \$$FREESURFER_HOME/bin/tcl_setup" >> $(DESTDIR)$(bindir)/unpacksdcmdir
	echo "unpacksdcmdir.tcl \$$argv" >> $(DESTDIR)$(bindir)/unpacksdcmdir
	chmod a+x $(DESTDIR)$(bindir)/unpacksdcmdir

EXTRA_DIST = $(tcl_SCRIPTS) tcl_setup vtk_setup kww_setup

# notes: 
# recon-all.2005 is the renamed recon-all-nmr, which is no longer used.
# recon-all is the official recon-all.
# recon-all-dev is the 'development' version of recon-all, which,
#               when changes are tested, get merged into recon-all.


# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
