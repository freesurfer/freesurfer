#
# Makefile.am
#

foodir=$(prefix)/average/Choi_JNeurophysiol12_MNI152

foo_datafiles=Choi2012_17Networks_ColorLUT.txt \
Choi2012_17NetworksConfidence_MNI152_FreeSurferConformed1mm_LooseMask.nii.gz \
Choi2012_17NetworksConfidence_MNI152_FreeSurferConformed1mm_TightMask.nii.gz \
Choi2012_17Networks_MNI152_FreeSurferConformed1mm_LooseMask.nii.gz \
Choi2012_17Networks_MNI152_FreeSurferConformed1mm_TightMask.nii.gz \
Choi2012_7Networks_ColorLUT.txt \
Choi2012_7NetworksConfidence_MNI152_FreeSurferConformed1mm_LooseMask.nii.gz \
Choi2012_7NetworksConfidence_MNI152_FreeSurferConformed1mm_TightMask.nii.gz \
Choi2012_7Networks_MNI152_FreeSurferConformed1mm_LooseMask.nii.gz \
Choi2012_7Networks_MNI152_FreeSurferConformed1mm_TightMask.nii.gz \
Choi_JNeurophysiol12_MNI152_README \
FSL_MNI152_FreeSurferConformed_1mm.nii.gz

install-data-local:
	test -z $(DESTDIR)$(foodir) || $(mkdir_p) $(DESTDIR)$(foodir)
	$(INSTALL_DATA) $(foo_datafiles) $(DESTDIR)$(foodir)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
