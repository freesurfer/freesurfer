# This makefile fragment runs the "check-TESTS" target from the generated
# Makefiles.  The idea is to allow developers to run make "check" targets
# outside of the nightly builds.  Tests have been empirically grouped 
# according to time needed to run into short, medium and long sets.
#
# To see the subdirectories with check-TESTS to run, and the test sets do,
# $ make -f check.mk
# or
# $ make -f check.mk list
#
# To run the check-TEST target for an indivudual directory, use the subdirectory name(s)
# as the target(s) to make,
# $ make -f check.mk mri_ca_register
# $ make -j8 -f check.mk mri_ca_register mrisp_write
#
# To run thre predefined groups/lists of tests from below do,
# $ make -f check.mk short	<--- runs short length tests serially, ~10-15 mins*
# $ make -f check.mk medium	<--- runs medium length tests in parallel with -j8, ~1 hour*
# $ make -f check.mk long	<--- runs long length tests in parallel with -j8, more than 1 hour*
#
# A good place to start is with,
#
# $ make -f check.mk short
#
# * machine configuration:
# I7 chip with 4 physcial 3.0 GHz cores (hyperthreading turned on
# to yield 8 virtual cores) with 16GB RAM using an SSD disk for I/O.

SHELL:=/bin/sh
MAKE:=/usr/bin/make

ifeq (,$(FREESURFER_HOME))
$(error Please set FREESURFER_HOME in the environment)
else
$(warning FREESURFER_HOME = $(FREESURFER_HOME))
endif

ifeq (,$(SUBJECTS_DIR))
$(error Please set SUBJECTS_DIR in the environment)
else
$(warning SUBJECTS_DIR = $(SUBJECTS_DIR))
endif

ifeq (,$(FS_LICENSE))
$(error Please set FS_LICENSE in the environment)
else
$(warning FS_LICENSE = $(FS_LICENSE))
endif

linux_default_subjects_dir:=
darwin_default_subjects_dir:=/Applications/freesurfer/subjects

os:=$(shell uname -s | tr -s '[A-Z]' '[a-z]')
ifeq ($(os),linux)
   ifeq (,$(SUBJECTS_DIR))
       SUBJECTS_DIR:=$(linux_default_subjects_dir)
   endif
   platform:=x86_64-$(os)
   platform_short:=linux
else ifeq ($(os),darwin)
   ifeq (,$(SUBJECTS_DIR))
       SUBJECTS_DIR:=$(darwin_default_subjects_dir)
   endif
   platform:=i386-mac
   platform_short:=mac
else
   $(error Cannot build for unknown os = $(os))
endif

ifeq (,$(SUBJECTS_DIR))
   $(error Please set SUBJECTS_DIR for $(platform_short) in environment prior to running tests)
endif
ifeq (,$(wildcard $(SUBJECTS_DIR)))
   $(error Cannot stat SUBJECTS_DIR=$(SUBJECTS_DIR) for $(platform_short))
endif
$(warning Using SUBJECTS_DIR=$(SUBJECTS_DIR))
export SUBJECTS_DIR

this_makefile:=$(lastword $(MAKEFILE_LIST))
this_makefile_path:=$(dir $(abspath $(this_makefile)))
# absolute path
# src_base:=$(this_makefile_path)
# relative path (know we are at the top of the tree)
src_base:=.
# $(warning src_base=$(src_base))

# get all Makefiles with generated check-TESTS targets
check_tests_dirs:=$(shell find $(src_base) -maxdepth 2 -name Makefile -exec grep "check-TESTS:" {} \; -print | grep Makefile | sed 's;^\.\/;;' | sed 's;\/\/;\/;' | sed 's;\/Makefile;;')
# $(warning check_tests_dirs=$(check_tests_dirs))

no_tests:=mris_extract_patches

# remove anything in no_tests
check_tests_filt1:=$(filter-out $(no_tests),$(check_tests_dirs))
# $(warning check_test_filt1=$(check_tests_filt1))

# manually ran tests to find those that take relatively longer
medium_tests:=mri_aparc2aseg mri_ca_normalize mri_ca_register mri_convert \
             mri_coreg mri_em_register mri_normalize mri_robust_register mris_anatomical_stats \
             mris_ca_label mris_convert mris_curvature mris_expand mris_fix_topology mris_inflate \
             mris_register mris_sphere mris_transform mris_volmask talairach_avi xml2
# $(warning medium=$(medium_tests))

long_tests = mri_ca_label mris_make_surfaces
# $(warning long=$(long_tests))

short_tests_filt1:=$(filter-out $(long_tests),$(check_tests_filt1))
short_tests:=$(filter-out $(medium_tests),$(short_tests_filt1))


.PHONY: list all_dirs
list:
	@echo "NO TESTS:"
	@echo  $(no_tests) | tr -s ' ' '\n'
	@echo
	@echo "LONG TESTS:"
	@echo  $(long_tests) | tr -s ' ' '\n'
	@echo
	@echo "MEDIUM TESTS:"
	@echo  $(medium_tests) | tr -s ' ' '\n'
	@echo
	@echo "SHORT TESTS:"
	@echo  $(short_tests) | tr -s ' ' '\n'
all_dirs:
	@echo "ALL check-TEST SUBDIRS:"
	@echo  $(check_tests_dirs) | tr -s ' ' '\n'
	@echo

# Setup as a target here each subdir with (generated) Makefile contaning check-TEST target.
# For cmake generated check targets, the list of subdirs $(check_tests_dirs) will only contain
# an entry for subdirs with a Makefile contaning that target. So if setup_configure and configure
# are not run to create Makefiles, then the empty list will cause make to exit w/o error (nothing to do).

.PHONY: $(check_tests_dirs)
$(check_tests_dirs):
	@echo "Making \"check-TESTS\" target under $@"
	(cd $@ && make && $(MAKE) check-TESTS)

# do not use := here as {short,medium,long}_tests may be appended to below
make_long_j8=$(MAKE) -k -j8 -f $(this_makefile) $(long_tests)
make_medium_j8=$(MAKE) -k -j8 -f $(this_makefile) $(medium_tests)
make_short_j1=$(MAKE) -k -j1 -f $(this_makefile) $(short_tests)

.PHONY: long medium short
.PHONY: long_tests medium_tests short_tests
long:
	@date
	$(make_long_j8)
	@date
medium:
	@date
	$(make_medium_j8)
	@date
short:
	@date
	$(make_short_j1)
	@date

# ==========================================================
# programmed tests for dirs w/o generated check-TESTS target

# FORCE:

define test_rule
.PHONY: $(1)
$(1):
	@echo
	@echo "############ Starting $(1) test ############"
	which $(2) && rm -rf $(3) && mkdir -p $(3)
ifeq ("False","$(5)")
	(cd $(1) && $(4))
else ifeq ("True","$(5)")
	@echo " >>>>> Test set to force exit status 0 <<<<< "
	(cd $(1) && $(4); exit 0)
endif
endef

define test_rule_no_cd
.PHONY: $(1)
$(1):
	@echo
	@echo "############ Starting $(1) test ############"
	which $(2) && rm -rf $(3) && mkdir -p $(3)
ifeq ("False","$(5)")
	$(4)
else ifeq ("True","$(5)")
	@echo " >>>>> Test set to force exit status 0 <<<<< "
	$(4); exit 0
endif
endef
# use $SUBJECTS_DIR data already in tree

local_test_base:=/tmp/test

# mri_binarize
#
# mri_binarize.bin --i ../mri_ca_label/testdata/mri/aseg.mgz --o /tmp/test/mri_binarize/mri_binarize.mgz --match 17
define test_mri_binarize
binarize_name := $1
binarize_prog := $$(binarize_name).bin
binarize_input := ../mri_ca_label/testdata/mri/aseg.mgz --match 17
binarize_output_dir := $(local_test_base)/$$(binarize_name)
binarize_output_file := $$(binarize_output_dir)/$$(binarize_name).mgz
binarize_cmd := $$(binarize_prog) --i $$(binarize_input) --o $$(binarize_output_dir)/$$(binarize_output_file)
binarize_force_exit_zero := False
$(call test_rule,$$(binarize_name),$$(binarize_prog),$$(binarize_output_dir),$$(binarize_cmd),$$(binarize_force_exit_zero))
endef
short_tests += mri_binarize

# mri_annotation2label
#
# mri_annotation2label --subject bert --hemi rh --labelbase /tmp/test/mri_annotation2label/labels/aparc-rh

define test_mri_annotation2label
anntl_name := $1
anntl_prog := $$(anntl_name)
anntl_input := --subject bert --hemi rh --labelbase
anntl_output_dir := $(local_test_base)/$$(anntl_name)/labels
anntl_output_file :=
anntl_cmd := $$(anntl_prog) $$(anntl_input) $$(anntl_output_dir)/aparc-rh
anntl_force_exit_zero := False
$(call test_rule,$$(anntl_name),$$(anntl_prog),$$(anntl_output_dir),$$(anntl_cmd),$$(anntl_force_exit_zero))
endef
short_tests += mri_annotation2label

# mri_concat
#
# mri_concat $SUBJECTS_DIR/fsaverage_sym/std.rh.*.mgh --o /tmp/rhout.mgh

define test_mri_concat
concat_name := $1
concat_prog := $$(concat_name)
concat_input := $(SUBJECTS_DIR)/fsaverage_sym/surf/std.rh.*.mgh
concat_output_dir := $(local_test_base)/$$(concat_name)
concat_output_file := rhout.mgh
concat_cmd := $$(concat_prog) $$(concat_input) --o $$(concat_output_dir)/$$(concat_output_file)
concat_force_exit_zero := False
$(call test_rule,$$(concat_name),$$(concat_prog),$$(concat_output_dir),$$(concat_cmd),$$(concat_force_exit_zero))
endef
short_tests += mri_concat

# mri_diff
#
# mri_diff ./fsaverage_sym/mri/lh.ribbon.mgz ./fsaverage_sym/mri/rh.ribbon.mgz --log /tmp/diff1

define test_mri_diff
diff_name := $1
diff_prog := $$(diff_name)
diff_input := $(SUBJECTS_DIR)/fsaverage_sym/mri/lh.ribbon.mgz $(SUBJECTS_DIR)/fsaverage_sym/mri/rh.ribbon.mgz
diff_output_dir := $(local_test_base)/$$(diff_name)
diff_output_file := diff_lh_and_rh_ribbon
diff_cmd := $$(diff_prog) $$(diff_input) --log $$(diff_output_dir)/$$(diff_output_file)
diff_force_exit_zero := True
$(call test_rule,$$(diff_name),$$(diff_prog),$$(diff_output_dir),$$(diff_cmd),$$(diff_force_exit_zero))
endef
short_tests += mri_diff

# mri_fuse_segmentations

# mri_info
#
# find $SUBJECTS_DIR -name "*.mgz" -exec mri_info {} \;

define test_mri_info
info_name := $1
info_prog := $$(info_name)
info_input := $(SUBJECTS_DIR) -name "*.mgz" -exec $$(info_name) {} \;
info_output_dir := $(local_test_base)/$$(info_name)
info_output_file := mri_info_output
info_cmd := find $$(info_input)
info_force_exit_zero := False
$(call test_rule,$$(info_name),$$(info_prog),$$(info_output_dir),$$(info_cmd),$$(info_force_exit_zero))
endef
short_tests += mri_info

# mri_label2label
#
# mri_label2label --srclabel $SUBJECTS_DIR/bert/label/lh.BA1_exvivo.label  --srcsubject bert --trglabel $SUBJECTS_DIR /cvs_avg35/label/lh.BA1.label --trgsubject cvs_avg35  --regmethod surface --hemi lh

define test_mri_label2label
l2l_name := $1
l2l_prog := $$(l2l_name)
l2l_input := --srclabel $(SUBJECTS_DIR)/bert/label/lh.BA1_exvivo.label  --srcsubject bert --trglabel $(SUBJECTS_DIR)/cvs_avg35/label/lh.BA1.label --trgsubject cvs_avg35  --regmethod surface --hemi lh
l2l_output_dir := $(local_test_base)/$$(l2l_name)
l2l_output_file := l2l_lh_and_rh_ribbon
l2l_cmd := $$(l2l_prog) $$(l2l_input)
l2l_force_exit_zero := False
$(call test_rule,$$(l2l_name),$$(l2l_prog),$$(l2l_output_dir),$$(l2l_cmd),$$(l2l_force_exit_zero))
endef
short_tests += mri_label2label

# mri label2vol
#
# mri_label2vol --label $SUBJECTS_DIR/bert/label/lh.cortex.label --regheader $SUBJECTS_DIR/bert/mri/rawavg.mgz --temp $SUBJECTS_DIR/bert/mri/aparc+aseg.mgz --o /tmp/l2v.nii

define test_mri_label2vol
l2v_name := $1
l2v_prog := $$(l2v_name)
l2v_input := --label $(SUBJECTS_DIR)/bert/label/lh.cortex.label --regheader $(SUBJECTS_DIR)/bert/mri/rawavg.mgz --temp $(SUBJECTS_DIR)/bert/mri/aparc+aseg.mgz
l2v_output_dir := $(local_test_base)/$$(l2v_name)
l2v_output_file := l2v.nii
l2v_cmd := $$(l2v_prog) $$(l2v_input) --o $$(l2v_output_dir)/$$(l2v_output_file)
l2v_force_exit_zero := False
$(call test_rule,$$(l2v_name),$$(l2v_prog),$$(l2v_output_dir),$$(l2v_cmd),$$(l2v_force_exit_zero))
endef
short_tests += mri_label2vol

# mri_log_likelihood

# mri_matrix_multiply
#
# mri_matrix_multiply -im $SUBJECTS_DIR/bert/mri/transforms/talairach.xfm -im $SUBJECTS_DIR/cvs_avg35/mri/transforms/talairach.xfm -v -om /tmp/mmult.xfm

define test_mri_matrix_multiply
mmult_name := $1
mmult_prog := $$(mmult_name)
mmult_input_1 := $(SUBJECTS_DIR)/bert/mri/transforms/talairach.xfm
mmult_input_2 := $(SUBJECTS_DIR)/cvs_avg35/mri/transforms/talairach.xfm
mmult_input := -v -im $$(mmult_input_1) -im $$(mmult_input_2)
mmult_output_dir := $(local_test_base)/$$(mmult_name)
mmult_output_file := mmult.xfm
mmult_result := $$(mmult_output_dir)/$$(mmult_output_file)
# output switch is -om !
mmult_cmd := cat $$(mmult_input_1) && cat $$(mmult_input_2) && $$(mmult_prog) $$(mmult_input) -om $$(mmult_result) && cat $$(mmult_result)
mmult_force_exit_zero := False
$(call test_rule,$$(mmult_name),$$(mmult_prog),$$(mmult_output_dir),$$(mmult_cmd),$$(mmult_force_exit_zero))
endef
short_tests += mri_matrix_multiply

# mri_normalize_tp2

# mri_nu_correct.mni
#
# mri_nu_correct.mni --no-rescale --i $SUBJECTS_DIR/bert/mri/orig.mgz --proto-iters 1000 --distance 50 --o /tmp/orig_nu.mgz

define test_mri_nu_correct.mni
mncm_name := $1
mncm_prog := $$(mncm_name)
mncm_input := --no-rescale --i $(SUBJECTS_DIR)/bert/mri/orig.mgz --proto-iters 1000 --distance 50  
mncm_output_dir := $(local_test_base)/$$(mncm_name)
mncm_output_file := orig_nu.mgz
mncm_result := $$(mncm_output_dir)/$$(mncm_output_file)
mncm_cmd := $$(mncm_prog) $$(mncm_input) --o $$(mncm_result)
mncm_force_exit_zero := False
# no such subdir as mri_nu_correct.mni
$(call test_rule_no_cd,$$(mncm_name),$$(mncm_prog),$$(mncm_output_dir),$$(mncm_cmd),$$(mncm_force_exit_zero))
endef
short_tests += mri_nu_correct.mni

# mri_pretess
#
# mri_pretess -test -w $SUBJECTS_DIR/bert/mri/wm.mgz wm $SUBJECTS_DIR/bert/mri/norm.mgz /tmp/wm_new.mgz

define test_mri_pretess
mrip_name := $1
mrip_prog := $$(mrip_name)
mrip_input := -test -w $(SUBJECTS_DIR)/bert/mri/wm.mgz wm $(SUBJECTS_DIR)/bert/mri/norm.mgz
mrip_output_dir := $(local_test_base)/$$(mrip_name)
mrip_output_file := wm_new.mgz
mrip_result := $$(mrip_output_dir)/$$(mrip_output_file)
mrip_cmd := $$(mrip_prog) $$(mrip_input) --o $$(mrip_result)
mrip_force_exit_zero := False
$(call test_rule_no_cd,$$(mrip_name),$$(mrip_prog),$$(mrip_output_dir),$$(mrip_cmd),$$(mrip_force_exit_zero))
endef
short_tests += mri_pretess

# mri_relabel_hypointensities
#
# mri_relabel_hypointensities $SUBJECTS_DIR/bert/mri/aseg.presurf.mgz $SUBJECTS_DIR/bert/surf /tmp/aseg_hypoint_out.mgz

define test_mri_relabel_hypointensities
mrrhypo_name := $1
mrrhypo_prog := $$(mrrhypo_name)
mrrhypo_input := $(SUBJECTS_DIR)/bert/mri/aseg.presurf.mgz $(SUBJECTS_DIR)/bert/surf
mrrhypo_output_dir := $(local_test_base)/$$(mrrhypo_name)
mrrhypo_output_file := aseg_hypoint_out.mgz
mrrhypo_result := $$(mrrhypo_output_dir)/$$(mrrhypo_output_file)
# no --o switch for output
mrrhypo_cmd := $$(mrrhypo_prog) $$(mrrhypo_input) $$(mrrhypo_result)
mrrhypo_force_exit_zero := False
$(call test_rule_no_cd,$$(mrrhypo_name),$$(mrrhypo_prog),$$(mrrhypo_output_dir),$$(mrrhypo_cmd),$$(mrrhypo_force_exit_zero))
endef
short_tests += mri_relabel_hypointensities

# mri_relabel_nonwm_hypos

# mri_robust_template
#
# Note addition of --debug generates output saying no permission to write (but not where), and exist status still zero
# Has exxample output for registration of TP1 to TP2
#
# mri_robust_template --debug -mov $SUBJECTS_DIR/bert/mri/orig/001.mgz $SUBJECTS_DIR/bert/mri/orig/002.mgz --average 1 --template /tmp/rawavg.mgz --satit --inittp 1 --fixtp --noit --iscale --subsample 200

define test_mri_robust_template
mrt_name := $1
mrt_prog := $$(mrt_name)
mrt_input := --debug -mov $(SUBJECTS_DIR)/bert/mri/orig/001.mgz $(SUBJECTS_DIR)/bert/mri/orig/002.mgz --average 1 --satit --inittp 1 --fixtp --noit --iscale --subsample 200
mrt_output_dir := $(local_test_base)/$$(mrt_name)
mrt_output_file := rawavg.mgz
mrt_result := $$(mrt_output_dir)/$$(mrt_output_file)
# no --o switch for output
mrt_cmd := $$(mrt_prog) $$(mrt_input) --template $$(mrt_result)
mrt_force_exit_zero := False
$(call test_rule_no_cd,$$(mrt_name),$$(mrt_prog),$$(mrt_output_dir),$$(mrt_cmd),$$(mrt_force_exit_zero))
endef
short_tests += mri_robust_template

# mri_seg_diff
#
# mri_seg_diff --debug --seg1 $SUBJECTS_DIR/bert/mri/aseg.auto.mgz --seg2 $SUBJECTS_DIR/bert/mri/aseg.mgz --diff /tmp/aseg_diff_out.mgz

define test_mri_seg_diff
msd_name := $1
msd_prog := $$(msd_name)
msd_input := --debug --seg1 $(SUBJECTS_DIR)/bert/mri/aseg.auto.mgz --seg2 $(SUBJECTS_DIR)/bert/mri/aseg.mgz
msd_output_dir := $(local_test_base)/$$(msd_name)
msd_output_file := aseg_diff_out.mgz
msd_result := $$(msd_output_dir)/$$(msd_output_file)
# no --o switch for output
msd_cmd := $$(msd_prog) $$(msd_input) --diff $$(msd_result)
msd_force_exit_zero := False
$(call test_rule_no_cd,$$(msd_name),$$(msd_prog),$$(msd_output_dir),$$(msd_cmd),$$(msd_force_exit_zero))
endef
short_tests += mri_seg_diff

# mri_segment
#
# mri_segment $SUBJECTS_DIR/bert/mri/brainmask.mgz /tmp/wm_segment_out.mgz -wlo 80

define test_mri_segment
mseg_name := $1
mseg_prog := $$(mseg_name)
mseg_input := $(SUBJECTS_DIR)/bert/mri/brainmask.mgz
mseg_output_dir := $(local_test_base)/$$(mseg_name)
mseg_output_file := wm_segment_out.mgz
mseg_result := $$(mseg_output_dir)/$$(mseg_output_file)
# no --o switch for output, cannot re-order -wlo 80 to be before output volume
mseg_cmd := $$(mseg_prog) $$(mseg_input) $$(mseg_result) -wlo 80
mseg_force_exit_zero := False
$(call test_rule_no_cd,$$(mseg_name),$$(mseg_prog),$$(mseg_output_dir),$$(mseg_cmd),$$(mseg_force_exit_zero))
endef
short_tests += mri_segment

# mri_stats2seg

# mri_surf2surf
# fails...
# mri_surf2surf --hemi lh --srcsubject bert --srcsurfval thickness --src_type curv --trgsubject bert --trgicoorder 7 --trgsurfval /tmp/bert-thickness-lh.img --trg_type analyze4d

define test_mri_surf2surf
ms2s_name := $1
ms2s_prog := $$(ms2s_name)
ms2s_input := --hemi lh --srcsubject bert --srcsurfval thickness --src_type curv --trgsubject bert --trgicoorder 7 --trg_type analyze4d
ms2s_output_dir := $(local_test_base)/$$(ms2s_name)
ms2s_output_file := bert-thickness-lh.img
ms2s_result := $$(ms2s_output_dir)/$$(ms2s_output_file)
# no --o switch for output
ms2s_cmd := $$(ms2s_prog) $$(ms2s_input) --trgsurfval $$(ms2s_result)
ms2s_force_exit_zero := False
$(call test_rule_no_cd,$$(ms2s_name),$$(ms2s_prog),$$(ms2s_output_dir),$$(ms2s_cmd),$$(ms2s_force_exit_zero))
endef
short_tests += mri_surf2surf

# mri_surf2vol

# mri_surfcluster

# mri_vol2surf

# mri_voltovol

# mri_voldiff

# mris_diff

# mris_divide_parcellation

# mris_label2annot

# mris_surface_stats

# mris_thickness

# mris_thickness_diff

# mris_topo_fixer

# pctsurfcon

# vertexvol

added_dirs = mri_binarize mri_annotation2label mri_concat mri_diff mri_info mri_label2label mri_label2vol mri_matrix_multiply mri_nu_correct.mni mri_pretess mri_relabel_hypointensities mri_robust_template mri_seg_diff mri_segment mri_surf2surf
# added_dirs = mri_surf2surf

.PHONY: $(added_dirs)
# create rules for test functions written above
$(foreach dir,$(added_dirs),$(eval $(call test_$(dir),$(dir))))
# for now, adding in tests with functions above to short tests (as don't take too long)
short_tests += $(added_dirs)
# $(warning short=$(short_tests))

