#
# Makefile.am
#

foodir=$(prefix)/average/Yeo_JNeurophysiol11_MNI152

foo_datafiles=Yeo_JNeurophysiol11_MNI152_README FSL_MNI152_FreeSurferConformed_1mm.nii.gz \
Yeo2011_17Networks_ColorLUT.txt Yeo2011_17Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii.gz \
Yeo2011_17Networks_MNI152_FreeSurferConformed1mm.nii.gz Yeo2011_7Networks_ColorLUT.txt \
Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii.gz Yeo2011_7Networks_MNI152_FreeSurferConformed1mm.nii.gz

install-data-local:
	test -z $(DESTDIR)$(foodir) || $(mkdir_p) $(DESTDIR)$(foodir)
	$(INSTALL_DATA) $(foo_datafiles) $(DESTDIR)$(foodir)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
