#
# Makefile.am for HippoSF (hippocampal subfield segmentation tools)
#

# platform-specific binaries found here:
SUBDIRS=linux_x86_64

# platform-independent data and scripts:
dist_bin_SCRIPTS = \
segmentSF_T1.sh \
segmentSF_T1T2.sh \
segmentSF_T2.sh

foodir=$(prefix)/average/HippoSF/atlas
foo_DATA=AtlasDump.mgz AtlasMesh.gz compressionLookupTable.txt

EXTRA_DIST=$(foo_DATA)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
