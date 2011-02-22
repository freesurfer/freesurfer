#
# Makefile.am
#

foodir=$(prefix)/docs
foo_DATA=\
	license.freesurfer_SLA.txt \
	license.minc.txt \
	license.fsl.txt \
	license.tcltk.txt \
	license.tix.txt \
	license.blt.txt \
	license.swig.txt \
	license.vxl.txt \
	license.thrust.txt \
	license.qt.html \
	license.vtk.txt \
	license.itk.txt \
	license.petsc.txt \
	license.boost.txt \
	license.nrrdio.txt \
	license.nifti.txt \
	license.gifti.txt \
	license.cuda.txt

EXTRA_DIST=$(foo_DATA)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
