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

if HAVE_OPENGL_LIBS
GLUT_THINGS=-I$(top_srcdir)/glut -DHAVE_GLUT_LIBS
else
GLUT_THINGS=
endif

if USE_LOCAL_TIFF
TIFF_THINGS=-I$(top_srcdir)/tiff
else
TIFF_THINGS=
endif

if USE_LOCAL_JPEG
JPEG_THINGS=-I$(top_srcdir)/jpeg
else
JPEG_THINGS=
endif

if USE_LOCAL_EXPAT
EXPAT_THINGS=-I$(top_srcdir)/expat
else
EXPAT_THINGS=
endif

NRRDIO_CPPFLAGS=-I$(top_srcdir)/include/NrrdIO \
	-DTEEM_DIO=$(TEEM_DIO) -DTEEM_32BIT=$(TEEM_32BIT) \
	-DTEEM_QNANHIBIT=$(TEEM_QNANHIBIT) -DTEEM_ENDIAN=$(TEEM_ENDIAN)

AM_CPPFLAGS=\
	-I$(top_srcdir)/include \
	-I$(top_srcdir)/include/dicom \
  -I$(top_srcdir)/xml2 \
	$(NRRDIO_CPPFLAGS) $(JPEG_THINGS) $(TIFF_THINGS) $(EXPAT_THINGS) \
	$(ITK_THINGS) $(NIFTI_THINGS) $(GLUT_THINGS) \
	$(GL_CFLAGS) $(VXL_CFLAGS) $(CUDA_CFLAGS) -DHAVE_ZLIB

AM_LDFLAGS=

LIBS=

noinst_LIBRARIES=libutils.a

libutils_main_SOURCES=\
	afni.c \
	annotation.c \
	art1.c \
	artmap.c \
	aseg_edit_svm.c \
	backprop.c \
	bfileio.c \
	box.c \
	Bruker.c \
	chklc.c \
	chronometer.c \
	cluster.c \
	class_array.c \
	cma.c \
	cmdargs.c \
	colortab.c \
	connectcomp.c \
	corio.c \
	cpputils.cpp \
	ctrpoints.c \
	cvector.c \
	dcdflib.c \
	dct.c \
	density.c \
	diag.c \
	DICOMRead.c \
	diffuse.c \
	dti.c \
	dtk.fs.c \
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
	fsPrintHelp.c \
	fsinit.c \
	gcaboundary.c \
	gca.c \
	gcamorph.c \
	gcarray.c \
	gclass.c \
	gcsa.c \
	getdelim.c \
	getline.c \
	gifti_io.c \
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
	ic42.c \
	ic12.c \
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
	MARS_DT_Boundary.c \
	matfile.c \
	matrix.c \
	mgh_filter.c \
	mgh_matrix.c \
	min_heap.c \
	mincutils.c \
	minmaxrc.c \
	morph.c \
	mosaic.c \
	mri2.c \
	mri.c \
  mriBSpline.c \
	mriclass.c \
	mri_conform.c \
	mricurv.c \
	mriflood.c \
	mriFunctionalDataAccess.c \
	mriHeadPointList.c \
	mrihisto.c \
	mri_identify.c \
	MRIio_old.c \
	mriio.c \
	mri_level_set.c \
	mrifilter.c \
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
	Progname.c \
	queue.c \
	randomfields.c \
	rbf.c \
	rec.c \
	region.c \
	registerio.c \
	resample.c \
	rescale.c \
	retinotopy.c \
	rfa.c \
	rforest.c \
	selxavgio.c \
	sig.c \
	signa.c \
	stats.c \
	stc.c \
	surfcluster.c \
	svm.c \
	tags.c \
	talairachex.c \
	thread.c \
	timer.c \
	transform.c \
	tritri.c \
	tukey.c \
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
	mriio_nrrd.c \
	mriio_nrrd_itk.cpp \
	numerics.cpp \
	chronometerpp.cpp\
	gcamorphtestutils.cpp\
	gcautils.cpp\
	gcalinearnode.cpp\
	gcalinearprior.cpp

libutils_cephes_SOURCES=\
	$(srcdir)/$(cephes_dir)/bdtr.c \
	$(srcdir)/$(cephes_dir)/chdtr.c \
	$(srcdir)/$(cephes_dir)/const.c \
	$(srcdir)/$(cephes_dir)/fdtr.c \
	$(srcdir)/$(cephes_dir)/gamma.c \
	$(srcdir)/$(cephes_dir)/igam.c \
	$(srcdir)/$(cephes_dir)/igami.c \
	$(srcdir)/$(cephes_dir)/incbet.c \
	$(srcdir)/$(cephes_dir)/incbi.c \
	$(srcdir)/$(cephes_dir)/mtherr.c \
	$(srcdir)/$(cephes_dir)/ndtr.c \
	$(srcdir)/$(cephes_dir)/ndtri.c \
	$(srcdir)/$(cephes_dir)/polevl.c \
	$(srcdir)/$(cephes_dir)/stdtr.c

libutils_nrrdio_SOURCES=\
	$(srcdir)/$(NrrdIO_dir)/754.c \
	$(srcdir)/$(NrrdIO_dir)/formatPNG.c \
	$(srcdir)/$(NrrdIO_dir)/accessors.c \
	$(srcdir)/$(NrrdIO_dir)/formatPNM.c \
	$(srcdir)/$(NrrdIO_dir)/array.c \
	$(srcdir)/$(NrrdIO_dir)/formatText.c \
	$(srcdir)/$(NrrdIO_dir)/arraysNrrd.c \
	$(srcdir)/$(NrrdIO_dir)/formatVTK.c \
	$(srcdir)/$(NrrdIO_dir)/axis.c \
	$(srcdir)/$(NrrdIO_dir)/gzio.c \
	$(srcdir)/$(NrrdIO_dir)/biffbiff.c \
	$(srcdir)/$(NrrdIO_dir)/keyvalue.c \
	$(srcdir)/$(NrrdIO_dir)/comment.c \
	$(srcdir)/$(NrrdIO_dir)/methodsNrrd.c \
	$(srcdir)/$(NrrdIO_dir)/defaultsNrrd.c \
	$(srcdir)/$(NrrdIO_dir)/miscAir.c \
	$(srcdir)/$(NrrdIO_dir)/dio.c \
	$(srcdir)/$(NrrdIO_dir)/mop.c \
	$(srcdir)/$(NrrdIO_dir)/encodingAscii.c \
	$(srcdir)/$(NrrdIO_dir)/parseAir.c \
	$(srcdir)/$(NrrdIO_dir)/encodingBzip2.c \
	$(srcdir)/$(NrrdIO_dir)/parseNrrd.c \
	$(srcdir)/$(NrrdIO_dir)/encoding.c \
	$(srcdir)/$(NrrdIO_dir)/encodingGzip.c \
	$(srcdir)/$(NrrdIO_dir)/encodingHex.c \
	$(srcdir)/$(NrrdIO_dir)/read.c \
	$(srcdir)/$(NrrdIO_dir)/encodingRaw.c \
	$(srcdir)/$(NrrdIO_dir)/reorder.c \
	$(srcdir)/$(NrrdIO_dir)/endianAir.c \
	$(srcdir)/$(NrrdIO_dir)/endianNrrd.c \
	$(srcdir)/$(NrrdIO_dir)/sane.c \
	$(srcdir)/$(NrrdIO_dir)/enum.c \
	$(srcdir)/$(NrrdIO_dir)/simple.c \
	$(srcdir)/$(NrrdIO_dir)/enumsNrrd.c \
	$(srcdir)/$(NrrdIO_dir)/string.c \
	$(srcdir)/$(NrrdIO_dir)/format.c \
	$(srcdir)/$(NrrdIO_dir)/subset.c \
	$(srcdir)/$(NrrdIO_dir)/formatEPS.c \
	$(srcdir)/$(NrrdIO_dir)/write.c \
	$(srcdir)/$(NrrdIO_dir)/formatNRRD.c

libutils_a_SOURCES=\
	$(libutils_main_SOURCES) \
	$(libutils_cephes_SOURCES) \
	$(libutils_nrrdio_SOURCES)

## ----
## CUDA
## ----

# BUILDCUDA is defined if configure.in finds CUDA
if BUILDCUDA
# rules for building cuda files
NVCCINCLUDES = $(CUDA_CFLAGS) -I$(top_srcdir)/include $(MNI_CFLAGS)

if HAVE_FERMI_GPU
# Fermi flags:
NVCC_ARCH = -DGCAMORPH_ON_GPU
CUDA_GCA_CFLAGS = -DGCAMORPH_ON_GPU
AM_CXXFLAGS = -DGCAMORPH_ON_GPU -DFS_CUDA -DCUDA_FORCE_SYNC -DCUDA_SHOW_TIMINGS
else
# Tesla flags (default):
NVCC_ARCH = 
CUDA_GCA_CFLAGS = 
endif

G++FLAGS = -Wall -Wextra -Weffc++

CUDA_DEPS = $(top_srcdir)/include/volumegpu.hpp\
            $(top_srcdir)/include/volumegpu_impl.hpp\
            $(top_srcdir)/include/mriframegpu.hpp\
            $(top_srcdir)/include/gcamorphgpu.hpp\
            $(top_srcdir)/include/mriconvolve_cuda.hpp\
            $(top_srcdir)/include/mrimean_cuda.hpp\
            $(top_srcdir)/include/gcamorphenergy.hpp\
            $(top_srcdir)/include/gcamorphtermgpu.hpp\
            $(top_srcdir)/include/fixedmap.hpp\
            $(top_srcdir)/include/vecvolgpu.hpp\
            $(top_srcdir)/include/ctfactory.hpp

.cu.o:
	$(NVCC) -o $@ -c $< $(NVCCFLAGS) $(NVCCINCLUDES) $(NVCC_ARCH) \
		-DFS_CUDA -DCUDA_FORCE_SYNC -DCUDA_SHOW_TIMINGS \
                --compiler-options -O3

noinst_LIBRARIES += libutils_cuda.a

libutils_cuda_a_SOURCES=\
	$(libutils_main_SOURCES) \
	$(libutils_cephes_SOURCES) \
	$(libutils_nrrdio_SOURCES) \
        mriframegpu.cu \
	devicemanagement.cu \
	mriconvolve_cuda.cu \
        mrivol2vol_cuda.cu \
	mrisurf_cuda.cu \
        affinegpu.cu\
        mrimean_cuda.cu\
        gcasgpu.cu\
        cudatypeutils.cu\
        gcamorphgpu.cu\
        gcamorphenergygpu.cu\
        gcamfots_cuda.cu\
        gcamorphtermgpu.cu\
        gcamorphcpu.cpp\
        gcamcomputegradient_cuda.cu\
        gcamregisterlevel_cuda.cu\
        gcamregisterpipeline_cuda.cu\
        mrilabels_cuda.cu\
        gcanodegpu.cu\
        gcapriorgpu.cu\
        vecvolgpu.cu

libutils_cuda_a_CFLAGS=-DFS_CUDA $(CUDA_CFLAGS) $(CUDA_GCA_CFLAGS)

# We have some templating dependencies
mriframegpu_cuda.o: $(CUDA_DEPS)
mriconvolve_cuda.o: $(CUDA_DEPS)
mrimean_cuda.o: $(CUDA_DEPS)
mrivol2vol_cuda.o: $(CUDA_DEPS)
affinegpu.o: $(CUDA_DEPS)
gcasgpu.o: $(CUDA_DEPS)
gcamorphgpu.o: $(CUDA_DEPS) $(top_srcdir)/include/volumegpucompare.hpp
cudatypeutils.o: $(top_srcdir)/include/cudatypeutils.hpp
gcamorphenergygpu.o: $(CUDA_DEPS)
gcamorphtermgpu.o: $(CUDA_DEPS)  $(top_srcdir)/include/volumecpu.hpp
gcamorphcpu.o: $(top_srcdir)/include/gcamorphcpu.hpp $(CUDA_DEPS)
mrilabels_cuda.o: $(CUDA_DEPS) $(top_srcdir)/include/mrilabels_cuda.hpp
gcanodegpu.o: $(CUDA_DEPS) $(top_srcdir)/include/gcanodegpu.hpp $(top_srcdir)/include/gcalinearnode.hpp
gcapriorgpu.o: $(CUDA_DEPS) $(top_srcdir)/include/gcapriorgpu.hpp $(top_srcdir)/include/gcalinearprior.hpp
gcamregisterlevel_cuda.o : $(CUDA_DEPS) $(top_srcdir)/include/gcamremovenegativenodes_cuda.hpp

CALLDEPTH ?= -1

%.calltree : %.dot
	dot -Grankdir=LR -Tpdf < $< > $*.pdf

%.dot : %.cflow
	cflow2dot < $< > $@

%.cflow : 
	cflow -d ${CALLDEPTH} \
              --format=posix --omit-arguments \
              --level-indent='0=\t' --level-indent='1=\t' \
              --level-indent=start='\t' --main=$* \
              *.c > $@


endif



bin_PROGRAMS=gifti_tool gifti_test fsPrintHelp xmlToHtml

gifti_tool_SOURCES=gifti_tool.c gifti_tool.h \
	gifti_io.c gifti_xml.c nifti1_io.c znzlib.c
gifti_tool_LDADD=$(top_builddir)/expat/libexpat.a -lz -lm $(LIBS_NIFTI)

gifti_test_SOURCES=gifti_test.c gifti_test.h gifti_io.c gifti_xml.c \
	nifti1_io.c znzlib.c
gifti_test_LDADD=$(top_builddir)/expat/libexpat.a -lz $(LIBS_NIFTI)

fsPrintHelp_SOURCES=fsPrintHelp.c
fsPrintHelp_CFLAGS=-DBUILD_MAIN
fsPrintHelp_LDADD=$(top_builddir)/xml2/libxml2.a -lz
fsPrintHelp_LDFLAGS=$(OS_LDFLAGS)

xmlToHtml_SOURCES=xmlToHtml.c
xmlToHtml_LDADD=$(top_builddir)/xml2/libxml2.a -lz

install-exec-hook:
	cp $(DESTDIR)$(bindir)/fsPrintHelp $(DESTDIR)$(bindir)/fshelp
	cp $(DESTDIR)$(bindir)/fsPrintHelp $(DESTDIR)$(bindir)/fsman

uninstall-hook:
	rm -f $(DESTDIR)$(bindir)/fs*

if USE_SSE_MATHFUN
check_PROGRAMS = sse_mathfun_test
sse_mathfun_test_SOURCES = sse_mathfun_test.c
sse_mathfun_test_LDADD = -lm

TESTS = sse_mathfun_test
endif


EXTRA_DIST=$(srcdir)/$(NrrdIO_dir)/copyright.NrrdIO $(CUDA_DEPS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
