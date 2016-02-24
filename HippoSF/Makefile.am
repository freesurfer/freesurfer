#
# Makefile.am for HippoSF (hippocampal subfield segmentation tools)
#

foodir=$(prefix)/average/HippoSF/atlas
foo_datafiles=AtlasDump.mgz AtlasMesh.gz compressionLookupTable.txt

if HAVE_MAC_OSX

SUBDIRS=mac_osx

else

SUBDIRS=linux_x86_64

endif

# platform-independent data and scripts:
dist_bin_SCRIPTS = \
segmentSF_T1.sh \
segmentSF_T1T2.sh \
segmentSF_T2.sh \
quantifyHippocampalSubfields.sh

install-data-local:
	test -z $(DESTDIR)$(foodir) || $(mkdir_p) $(DESTDIR)$(foodir)
	$(INSTALL_DATA) $(foo_datafiles) $(DESTDIR)$(foodir)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
