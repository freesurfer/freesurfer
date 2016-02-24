#
# Makefile.am for HippoSF/linux_x86_64
#

dist_bin_SCRIPTS = \
	run_segmentSubjectT1_autoEstimateAlveusML.sh \
	run_segmentSubjectT1T2_autoEstimateAlveusML.sh \
	run_segmentSubjectT2_autoEstimateAlveusML.sh

#
# We handle the  binary files differently because binary 
# files are treated as a separate compenent to the 
# repository, and thus may not be present at compile time.
#

foodir=$(prefix)/bin

foo_binaryfiles = \
	segmentSubjectT1_autoEstimateAlveusML \
	segmentSubjectT1T2_autoEstimateAlveusML \
	segmentSubjectT2_autoEstimateAlveusML

install-data-local:
	test -z $(DESTDIR)$(foodir) || $(mkdir_p) $(DESTDIR)$(foodir)
	$(INSTALL) -m 755 $(foo_binaryfiles) $(DESTDIR)$(foodir)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra

