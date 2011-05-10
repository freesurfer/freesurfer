#
# Makefile.am for scripts
#

dist_bin_SCRIPTS = \
kvlSegmentHippocampalSubfields.sh \
kvlQuantifyHippocampalSubfieldSegmentations.sh \
kvlTestHippocampalSubfields.sh

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
