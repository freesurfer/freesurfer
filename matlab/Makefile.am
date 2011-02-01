#
# Makefile.am for matlab files
#

mfilesdir=$(prefix)/matlab
siemensdir=$(prefix)/matlab

mfiles_DATA=\
	angles2rotmat.m \
	barsegstats.m \
	bmmcost.m \
	bmm_mcvect.m \
	bmm_mcvhist.m \
	bmmroc.m \
	cc_cut_afd.m \
	cc_cut_dir_afd.m \
	cc_cut_table.m \
	combine_labels.m \
	convert_fieldsign.m \
	convert_unwarp_resample.m \
	cortical_labeling_afd_txt.m \
	cortical_labeling_dir_afd_txt.m \
	cortical_labeling_table.m \
	dice_labels.m \
	dtifa.m \
	dtimatrix.m \
	fisher_twoclass.m \
	fread3.m \
	fwrite3.m \
	isdicomfile.m \
	juelichmat2mat.m \
	labelic.m \
	llbmm.m \
	load_analyze_hdr.m \
	load_analyze.m \
	load_cor.m \
	load_csd.m \
	load_dicom_fl.m \
	load_dicom_series.m \
	load_gca.m \
	load_ima.m \
	load_mgh2.m \
	load_mgh.m \
	load_nifti_hdr.m \
	load_nifti.m \
	load_segstats.m \
	lta_read.m \
	mksubfov.m \
	mri_cdf2p.m \
	MRIeuler2Mdc.m \
	MRIextractImage.m \
	MRIfspec.m \
	MRIisBHDR.m \
	MRIisMGH.m \
	mri_kurtosis.m \
	MRImdc2euler.m \
	MRIread.m \
	MRIseg2labelxyz.m \
	MRIsegReg.m \
	MRIvol2vol.m \
	MRIvote.m \
	MRIwrite.m \
	mri_zcdf.m \
	mris_display.m \
	peakfinder.m \
	pons_cut_afd.m \
	pons_cut_dir_afd.m \
	pons_cut_table.m \
	randb.m \
	read_all.m \
	read_annotation.m \
	read_ascii_curv.m \
	read_asc.m \
	read_cor.m \
	read_csf_patch.m \
	read_curv.m \
	read_eccen_patch.m \
	read_freq_patch.m \
	read_fscolorlut.m \
	read_genesis_image.m \
	read_label.m \
	read_label_old.m \
	read_moviebyu.m \
	read_patch.m \
	readrec.m \
	ReadSiemensPhysio.m\
	read_siemens_header.m \
	read_siemens_image.m \
	read_smooth_eccen.m \
	read_surf.m \
	read_type.m \
	read_vf.m \
	read_wfile.m \
	remove_spaces.m \
	ribbon_afd.m \
	ribbon_dir_afd.m \
	ribbon_table.m \
	rotmat.m \
	rotmat2angles.m \
	sampleSize.m \
	save_cor.m \
	save_mgh2.m \
	save_mgh.m \
	save_nifti.m \
	ssblochgrad.m \
	ssbloch.m \
	stringunique.m \
	strlen.m \
	subcortical_labeling_afd.m \
	subcortical_labeling_dir_afd.m \
	subcortical_labeling_table.m \
	surf_registration_afd.m \
	surf_registration_stats.m \
	surf_registration_table.m \
	talairaching_afd.m \
	talairaching_dir_afd.m \
	talairaching_stats.m \
	talairaching_table.m \
	unwarp_init_globals.m \
	unwarp_resample.m \
	unwarp_scanners_table.m \
	vox2ras_0to1.m \
	vox2ras_1to0.m \
	vox2ras_dfmeas.m \
	vox2ras_ksolve.m \
	vox2ras_rsolveAA.m \
	vox2ras_rsolve.m \
	vox2ras_tkreg.m \
	vox2rasToQform.m \
	wm_seg_afd.m \
	wm_seg_dir_afd.m \
	wm_seg_table.m \
	write_analyze_hdr.m \
	write_annotation.m \
	write_ascii_curv.m \
	write_curv.m \
	write_label.m \
	write_wfile.m \
	xfm_read.m

if NMR_INSTALL
siemens_DATA=\
	read_meas_dat.m \
	read_meas_prot.m \
	read_meas_prot__struct.m
else
siemens_DATA=
endif

EXTRA_DIST=$(mfiles_DATA)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=$(siemens_DATA)
include $(top_srcdir)/Makefile.extra
