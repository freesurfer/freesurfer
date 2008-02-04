#
# Makefile.am for libutils.a
#

NrrdIO_dir = NrrdIO
cephes_dir = cephes
test_dir = test

SUBDIRS=$(NrrdIO_dir) $(cephes_dir) $(test_dir)

if HAVE_ITK_LIBS
ITK_THINGS=$(ITK_CFLAGS) -DHAVE_ITK_LIBS
else
ITK_THINGS=
endif

if HAVE_NIFTI_LIBS
NIFTI_THINGS=$(NIFTI_CFLAGS) -DHAVE_NIFTI_LIBS
else
NIFTI_THINGS=
endif

if HAVE_GLUT_LIBS
GLUT_THINGS=$(GLUT_CFLAGS) -DHAVE_GLUT_LIBS
else
GLUT_THINGS=
endif

AM_CPPFLAGS=$(X_CFLAGS) $(VXL_CFLAGS) \
	-I$(top_srcdir)/include \
	-I$(top_srcdir)/include/dicom \
	-I$(top_srcdir)/include/NrrdIO \
	$(ITK_THINGS) $(NIFTI_THINGS) $(GLUT_THINGS) -DHAVE_ZLIB

AM_LDFLAGS=

noinst_LIBRARIES= libutils_main.a libutils.a

libutils_main_a_SOURCES=\
	afni.c \
	annotation.c \
	art1.c \
	artmap.c \
	backprop.c \
	bfileio.c \
	box.c \
	Bruker.c \
	chklc.c \
	cluster.c \
	cma.c \
	cmdargs.c \
	colortab.c \
	corio.c \
	ctrpoints.c \
	cvector.c \
	dcdflib.c \
	dct.c \
	density.c \
	diag.c \
	DICOMRead.c \
	diffuse.c \
	dti.c \
	eigen.c \
	error.c \
	evschutils.c \
	fftutils.c \
	field_code.c \
	filter.c \
	fio.c \
	flash.c \
	fmarchmesh.c \
	fmriutils.c \
	fsenv.c \
	fsglm.c \
	gcaboundary.c \
	gca.c \
	gcamorph.c \
	gcarray.c \
	gclass.c \
	gcsa.c \
	getdelim.c \
	getline.c \
	gifti.c \
	gifti_local.c \
	gifti_xml.c \
	gw_ic2562.c \
	gw_utils.c \
	handle.c \
	heap.c \
	hippo.c \
	histo.c \
	ic163842.c \
	ic2562.c \
	ic642.c \
	icosahedron.c \
	image.c \
	imageio.c \
	imautils.c \
	intensity_eig.c \
	inverse.c \
	ipmpar.c \
	kernel.c \
	kinput.c \
	label.c \
	lpafile.c \
	machine.c \
	matfile.c \
	matrix.c \
	mgh_filter.c \
	mgh_matrix.c \
	mincutils.c \
	minmaxrc.c \
	morph.c \
	mosaic.c \
	mri2.c \
	mri.c \
	mriclass.c \
	mri_conform.c \
	mricurv.c \
	mrifilter.c \
	mriflood.c \
	mriFunctionalDataAccess.c \
	mriHeadPointList.c \
	mrihisto.c \
	mri_identify.c \
	mriio.c \
	MRIio_old.c \
	mri_level_set.c \
	mrimorph.c \
	mrinorm.c \
	mripolv.c \
	mriprob.c \
	mrisbiorthogonalwavelets.c \
	mrisegment.c \
	mriset.c \
	mrishash.c \
	mrisp.c \
	mriSurface.c \
	mrisurf.c \
	mrisutils.c \
	mri_tess.c \
	mri_topology.c \
	mri_transform.c \
	mriTransform.c \
	mriVolume.c \
	mrivoxel.c \
	nifti1_io.c \
	offset.c \
	path.c \
	pdf.c \
	pgmstubs.c \
	prime.c \
	qdecutils.c \
	queue.c \
	randomfields.c \
	rbf.c \
	rec.c \
	region.c \
	registerio.c \
	resample.c \
	rescale.c \
	retinotopy.c \
	selxavgio.c \
	sig.c \
	signa.c \
	stats.c \
	stc.c \
	stub.c \
	surfcluster.c \
	svm.c \
	tags.c \
	talairachex.c \
	thread.c \
	timer.c \
	transform.c \
	tritri.c \
	tukey.c \
	unwarpGradientNonlinearity.c \
	utils.c \
	version.c \
	vlabels.c \
	volcluster.c \
	voxlist.c \
	x3DList.c \
	xDebug.c \
	xGLutWindow.c \
	xGrowableArray.c \
	xList.c \
	xSparseVolume.c \
	xTypes.c \
	xUndoList.c \
	xUtilities.c \
	xVoxel.c \
	znzlib.c \
	fs_cost_function.cpp \
	fs_lbfgs.cpp \
	fs_lbfgs_observer.cpp \
	fs_lbfgs_subject.cpp \
	fs_powell.cpp \
	mri_fastmarching.cpp \
	mriio_nrrd.cpp \
	numerics.cpp

libutils_a_SOURCES=

am_libutils_a_OBJECTS= \
	$(libutils_main_a_OBJECTS) \
	$(top_builddir)/utils/$(cephes_dir)/libcephes.a \
	$(top_builddir)/utils/$(NrrdIO_dir)/libnrrdio.a

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
