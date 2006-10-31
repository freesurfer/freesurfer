#
# Makefile.am
#
# license files should be copied over to the install directory

foodir=$(prefix)/docs
foo_DATA=license.freesurfer \
	license.freesurfer_opensource \
	license.minc \
	license.fsl \
	license.tcltk \
	license.tix \
	license.slicer \
	license.swig \
	license.vxl

EXTRA_DIST=$(foo_DATA)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
