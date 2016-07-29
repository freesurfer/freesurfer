
#
# Makefile.am for HippoSF/mac_osx
#

dist_bin_SCRIPTS = \
	run_segmentSubjectT1_autoEstimateAlveusML.sh \
	run_segmentSubjectT1T2_autoEstimateAlveusML.sh \
	run_segmentSubjectT2_autoEstimateAlveusML.sh \
	run_SegmentSubfieldsT1Longitudinal.sh

#
# We handle the  binary files differently because binary 
# files are treated as a separate compenent to the 
# repository, and thus may not be present at compile time.
#

foodir=$(prefix)/bin

install-data-local:
	cp -RL segmentSubjectT1_autoEstimateAlveusML.app $(foodir)
	chmod -R 755 $(foodir)/segmentSubjectT1_autoEstimateAlveusML.app
	cp -RL segmentSubjectT1T2_autoEstimateAlveusML.app $(foodir)
	chmod -R 755 $(foodir)/segmentSubjectT1T2_autoEstimateAlveusML.app
	cp -RL segmentSubjectT2_autoEstimateAlveusML.app $(foodir)
	chmod -R 755 $(foodir)/segmentSubjectT2_autoEstimateAlveusML.app
	cp -RL SegmentSubfieldsT1Longitudinal.app $(foodir)
	chmod -R 755 $(foodir)/SegmentSubfieldsT1Longitudinal.app

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
