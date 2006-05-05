#
# Makefile.am
#
# licence files should be copied over to the install directory
foodir=$(prefix)/docs
foo_DATA=license.freesurfer \
	citation.minc \
	copywrite.minc \
	license.minc \
	registration.notes \
	license.mesa \
	license.fsl \
	license.gnu \
	license.tcltk \
	license.html_lib \
	license.tix \
	license.afni \
	license.slicer \
	license.swig

EXTRA_DIST=$(foo_DATA)


# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
