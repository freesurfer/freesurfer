#
# Makefile.am
#

foodir=$(prefix)/average/Buckner_JNeurophysiol11_MNI152

foo_datafiles=Buckner2011_17Networks_ColorLUT.txt \
Buckner2011_17NetworksConfidence_MNI152_FreeSurferConformed1mm_LooseMask.nii.gz \
Buckner2011_17NetworksConfidence_MNI152_FreeSurferConformed1mm_TightMask.nii.gz \
Buckner2011_17Networks_MNI152_FreeSurferConformed1mm_LooseMask.nii.gz \
Buckner2011_17Networks_MNI152_FreeSurferConformed1mm_TightMask.nii.gz \
Buckner2011_7Networks_ColorLUT.txt \
Buckner2011_7NetworksConfidence_MNI152_FreeSurferConformed1mm_LooseMask.nii.gz \
Buckner2011_7NetworksConfidence_MNI152_FreeSurferConformed1mm_TightMask.nii.gz \
Buckner2011_7Networks_MNI152_FreeSurferConformed1mm_LooseMask.nii.gz \
Buckner2011_7Networks_MNI152_FreeSurferConformed1mm_TightMask.nii.gz \
FSL_MNI152_FreeSurferConformed_1mm.nii.gz \
Buckner_JNeurophysiol11_MNI152_README


install-data-local:
	test -z $(DESTDIR)$(foodir) || $(mkdir_p) $(DESTDIR)$(foodir)
	$(INSTALL_DATA) $(foo_datafiles) $(DESTDIR)$(foodir)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
