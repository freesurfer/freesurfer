#
# Makefile.am for scripts
#
dist_bin_SCRIPTS = \
archive_pull \
browse-minc-header.tcl \
browse-sessions \
check_subject \
cor_to_minc \
csurf \
findsession \
fix_subject \
fix_subject_corrected \
fix_subject_corrected-lh \
fix_subject_corrected-rh \
fix_subject-lh \
fix_subject_on_seychelles \
fix_subject-rh \
fsl_rigid_register \
grad_unwarp \
inflate_subject \
inflate_subject3 \
inflate_subject-lh \
inflate_subject_new \
inflate_subject_new-lh \
inflate_subject_new-rh \
inflate_subject-rh \
inflate_subject_sc \
label_child \
label_elderly_subject \
label_subject \
label_subject_flash \
label_subject_mixed \
make_final_surfaces_corrected_subject \
make_final_surfaces_subject \
make_final_surfaces_subject-lh \
make_final_surfaces_subject_on_seychelles \
make_final_surfaces_subject-rh \
make_final_surfaces_subject_richard \
make_quick_surfaces_subject \
make-segvol-table \
map_all_labels \
map_all_labels-lh \
map_central_sulcus \
minc2seqinfo \
mkfsdist \
mkima_index.tcl \
mkmnc_index.tcl \
mksubjdirs \
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
mri_autorecon \
mri-func2sph \
mri-funcvits \
mri_mergelabels \
mri_motion_correct \
mri_motion_correct2 \
mri-sph2surf \
mri-structvits \
ms_refine_subject \
parcellate_subject \
parcellate_subject-lh \
parcellate_subject-rh \
printdat \
progressbar.tcl \
recon-all-nmr \
recon-all \
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
segment_monkey \
segment_subject \
segment_subject_notal \
segment_subject_notal2 \
segment_subject_old_skull_strip \
segment_subject_sc \
segment_subject_talmgh \
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
check_func_ras_coords.tcl \
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


recon-all: recon-all-nmr
	cat ./recon-all-nmr | grep -v dASEGd > ./recon-all
	chmod a+x ./recon-all

# use , as the separator so that the directory path / can be passed
SEDFLAGS = -e 's,xEXEC_PREFIXx,$(exec_prefix),g' \
	   -e 's,xPREFIXx,$(prefix),g' \
	   -e 's,xTIXWISHx,$(TIXWISH),g'

csurf: csurf.in
	sed < $< $(SEDFLAGS) > $@ 
	chmod +x $@


SEDFLAGS2 = -e 's,xFREESURFER_HOMEx,$(prefix),g' \
	    -e 's,xMNI_BINDIRx,$(mni_bindir),g' \
	    -e 's,xTCL_LIBDIRx,$(tcl_libdir),g'

freesurferenv.csh: freesurferenv.xxx
	sed < $< $(SEDFLAGS2) > $@
	chmod +x $@

foodir=$(prefix)
foo_SCRIPTS=freesurferenv.csh

# work on tixwish files
fsgdfPlot.tcl : fsgdfPlot.xxx
	sed < $< -e 's,/usr/bin/tixwish,$(TIXWISH),' > $@

tkmedit.tcl : tkmedit.xxx
	sed < $< -e 's,/usr/bin/tixwish,$(TIXWISH),' > $@

tkm_functional.tcl: tkm_functional.xxx
	sed < $< -e 's,/usr/bin/tixwish,$(TIXWISH),' > $@

tkm_wrappers.tcl: tkm_wrappers.xxx
	sed < $< -e 's,/usr/bin/tixwish,$(TIXWISH),' > $@

tksurfer.tcl: tksurfer.xxx
	sed < $< -e 's,/usr/bin/tixwish,$(TIXWISH),' > $@
