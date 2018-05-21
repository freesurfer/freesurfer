#
# Makefile.am for scripts/InfantFS
#

dist_bin_SCRIPTS = \
atlassubjects.csh \
create_pial_surfaces_mprage_subject.csh \
create_wm_surfaces_mprage_subject.csh \
find_cerebellar_wm_labels.csh \
find_cerebral_gm_labels.csh \
find_cerebral_wm_labels.csh \
find_neighbors_byAge.csh \
find_neighbors_byAge_extended.csh \
find_neighbors_byMI.csh \
find_neighbors_byMI_extended.csh \
find_putamen_labels.csh \
infant_recon_all \
LabelFusion \
run_LabelFusion.sh \
run_LabelFusion_wrapper.csh \
set_babydev_packages.csh \
set_infant_recon_params.csh \
SetLabelFusionParams.csh\
subjects_example.csh

xmldir=$(prefix)/docs/xml
xml_DATA=

TESTS = 

EXTRA_DIST = 

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
