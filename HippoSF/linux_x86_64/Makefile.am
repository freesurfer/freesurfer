#
# Makefile.am for HippoSF
#

if HAVE_MAC_OSX
else

# these the runtime-compiled binaries (and wrappers) for Linux x86_64

dist_bin_SCRIPTS = \
segmentSubjectT1_autoEstimateAlveusML \
segmentSubjectT1T2_autoEstimateAlveusML \
segmentSubjectT2_autoEstimateAlveusML \
run_segmentSubjectT1_autoEstimateAlveusML.sh \
run_segmentSubjectT1T2_autoEstimateAlveusML.sh \
run_segmentSubjectT2_autoEstimateAlveusML.sh

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra

endif # install if not a Mac (default is Linux 64b)
