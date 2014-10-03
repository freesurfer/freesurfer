#
# Makefile.am
#

foodir=$(prefix)/data/GEMS

foo_datafiles=\
CurrentMeshCollection30.gz \
compressionLookupTable.txt \
compressionLookupTable_left.txt \
imageDump.mgz \
nonPartialVolumeStatsLeftReferenceForBert.txt \
nonPartialVolumeStatsRightReferenceForBert.txt \
kvlQuantifyHippocampalSubfieldSegmentations.sh.help.xml \
kvlSegmentHippocampalSubfields.sh.help.xml

install-data-local:
	test -z $(foodir) || $(mkdir_p) $(foodir)
	$(INSTALL_DATA) $(foo_datafiles) $(foodir)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
