#
# Makefile.am for BrainstemSS (brainstem substructure segmentation tools)
#

foodir=$(prefix)/average/BrainstemSS/atlas
foo_datafiles=AtlasDump.mgz AtlasMesh.gz compressionLookupTable.txt

if HAVE_MAC_OSX

SUBDIRS=mac_osx

else

SUBDIRS=linux_x86_64

endif

# platform-independent data and scripts:
dist_bin_SCRIPTS = \
segmentBS.sh \
quantifyBrainstemStructures.sh

install-data-local:
	test -z $(DESTDIR)$(foodir) || $(mkdir_p) $(DESTDIR)$(foodir)
	$(INSTALL_DATA) $(foo_datafiles) $(DESTDIR)$(foodir)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
