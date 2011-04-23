#
# Makefile.am for scripts
#


SUBDIRS=WholeBrainParcellation \
		exVivoHippocampi 

dist_bin_SCRIPTS = \
kvlcoregisterIntensityBased \
kvlnonPartialVolumeStatsRightReferenceForBert.txt \
kvlTestHippocampalSubfields.sh \
kvlcoregisterASEGBased \
kvltestLinearClassifier.m \
kvlanalyzeADNIResults.m \
kvlSegmentHippocampalSubfields.sh.help.txt \
kvldoAllCombinationsOfSpecialSubject.sh \
kvlanalyzeResults.m \
kvlinspectAllResults.sh \
kvldoSpecialSubject.sh \
kvlQuantifyHippocampalSubfieldSegmentations.sh.help.txt \
kvlkillAll.sh \
kvlpermutationTest.m \
kvlQuantifyHippocampalSubfieldSegmentations.sh \
kvlSegmentHippocampalSubfields.sh \
kvlnonPartialVolumeStatsLeftReferenceForBert.txt \
kvldoSubjectFSPipelineForNifti.sh \
kvlcoregisterFreeSurferLabels.sh \
kvlmakeVolumePlot.m \
kvldoAllSubjects.sh \
kvldoAllSubjectsFSPipelineForNifti.sh \
kvlcollectResults.sh

# ensure kvlSetHardCodedLocations.sh correctly points to stuff
install-exec-hook:
	echo "#!/bin/bash" > $(DESTDIR)$(bindir)/kvlSetHardCodedLocations.sh
	echo "PATH=\$$PATH:\$$FREESURFER_HOME/bin " >> $(DESTDIR)$(bindir)/kvlSetHardCodedLocations.sh
	echo "atlasDirectory=\$$FREESURFER_HOME/data/Atlas3D" >> $(DESTDIR)$(bindir)/kvlSetHardCodedLocations.sh
	chmod a+x $(DESTDIR)$(bindir)/kvlSetHardCodedLocations.sh

uninstall-hook:
	rm -f $(DESTDIR)$(bindir)/kvlSetHardCodedLocations.sh

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
