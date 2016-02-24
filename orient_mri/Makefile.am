
if USE_LOCAL_TIFF
TIFF_THINGS=-I$(top_srcdir)/tiff -L$(top_srcdir)/tiff
else
TIFF_THINGS=$(TIFF_CFLAGS)
endif

if USE_LOCAL_JPEG
JPEG_THINGS=-I$(top_srcdir)/jpeg -L$(top_srcdir)/jpeg
else
JPEG_THINGS=$(JPEG_CFLAGS)
endif


AM_CFLAGS=-I$(top_srcdir)/include \
	$(TIFF_THINGS) $(JPEG_THINGS) \
	$(VTK_CFLAGS) $(KWWIDGETS_CFLAGS) $(TCL_CFLAGS)
AM_CXXFLAGS=-Wno-deprecated -I$(top_srcdir)/include \
	$(TIFF_THINGS) $(JPEG_THINGS)\
	-I$(top_srcdir)/vtkfsio \
	-I$(top_srcdir)/vtkutils \
	$(VTK_CXXFLAGS) $(KWWIDGETS_CXXFLAGS) $(TCL_CFLAGS)

if ENABLE_TCLTK_APPS
bin_PROGRAMS = orient_mri
orient_mri_SOURCES= \
	OrientMRIMain.cxx \
	vtkKWOrientMRIApp.cxx vtkKWOrientMRIApp.h vtkKWOrientMRIAppTcl.cxx \
	vtkKWOrientMRIWindow.cxx vtkKWOrientMRIWindow.h \
	vtkKWOrientMRIWindowTcl.cxx \
	vtkKWOrientMRIView2D.cxx vtkKWOrientMRIView2D.h \
	vtkKWOrientMRIView2DTcl.cxx \
	vtkKWOrientMRIView3D.cxx vtkKWOrientMRIView3D.h \
	vtkKWOrientMRIView3DTcl.cxx \
	vtkOrientMRIInteractorStyleView3D.cxx \
	vtkOrientMRIInteractorStyleView3D.h \
	vtkOrientMRIInteractorStyleView2D.cxx \
	vtkOrientMRIInteractorStyleView2D.h \
	OrientMRILibInit.cxx
orient_mri_LDADD= \
	$(top_builddir)/vtkfsio/libvtkfsio.a \
	$(top_builddir)/vtkutils/libvtkutils.a \
	$(addprefix $(top_builddir)/, $(LIBS_MGH)) $(LIBS_TCL_OPENGL) \
	$(VTK_LIBS) $(VTK_VERDICT_LIB) -lvtkWidgets -lvtkGraphics \
	-lvtkImaging -lvtkFiltering \
	-lvtkCommon -lvtksys \
	-lvtkHybrid $(VTK_HDF5_LIB) \
	-lvtkexoIIc -lvtkNetCDF $(VTK_NETCDFCXX_LIB) \
	-lvtkVolumeRendering -lvtkRendering -lvtkftgl \
	$(VTK_EXPAT_LIB) $(VTK_FREETYPE_LIB) \
	$(VTK_TIFF_LIB) $(VTK_JPEG_LIB) $(VTK_SQLITE_LIB) \
	$(VTK_PNG_LIB) $(VTK_Z_LIB) $(VTK_METAIO_LIB) \
	-lvtkIO	$(VTK_LSDYNA_LIB) -lvtkDICOMParser $(VTK_MPEG2ENCODE_LIB) \
	-lvtkImagingTCL -lvtkVolumeRenderingTCL -lvtkRenderingTCL \
	-lvtkFilteringTCL -lvtkWidgetsTCL \
	-lvtkHybridTCL -lvtkGraphicsTCL -lvtkImagingTCL -lvtkIOTCL \
	-lvtkCommonTCL \
	$(KWWIDGETS_LIBS) -lKWWidgets
orient_mri_LDFLAGS=$(OS_LDFLAGS) $(LDFLAGS_TCL)

orient_mridir=$(prefix)/lib/resource
orient_mri_DATA=OrientMRIIcons.txt

vtkKWOrientMRIAppTcl.cxx: vtkKWOrientMRIApp.h
	${VTK_DIR}/bin/vtkWrapTcl vtkKWOrientMRIApp.h \
	${VTK_LIB_VTK}/hints 1 vtkKWOrientMRIAppTcl.cxx

vtkKWOrientMRIWindowTcl.cxx: vtkKWOrientMRIWindow.h
	${VTK_DIR}/bin/vtkWrapTcl vtkKWOrientMRIWindow.h \
	${VTK_LIB_VTK}/hints 1 vtkKWOrientMRIWindowTcl.cxx

vtkKWOrientMRIView2DTcl.cxx: vtkKWOrientMRIView2D.h
	${VTK_DIR}/bin/vtkWrapTcl vtkKWOrientMRIView2D.h \
	${VTK_LIB_VTK}/hints 1 vtkKWOrientMRIView2DTcl.cxx

vtkKWOrientMRIView3DTcl.cxx: vtkKWOrientMRIView3D.h
	${VTK_DIR}/bin/vtkWrapTcl vtkKWOrientMRIView3D.h \
	${VTK_LIB_VTK}/hints 1 vtkKWOrientMRIView3DTcl.cxx

OrientMRILibInit.cxx: OrientMRILibInit.data
	${VTK_DIR}/bin/vtkWrapTclInit OrientMRILibInit.data OrientMRILibInit.cxx

# put a wrapper around the bin, used to setup tcltktixblt,vtk,kww enviro vars
install-exec-hook:orient_mri
	cp $(top_builddir)/orient_mri/orient_mri $(DESTDIR)$(bindir)/orient_mri.bin
	echo "#!/bin/tcsh -ef" > $(DESTDIR)$(bindir)/orient_mri
	echo "source \$$FREESURFER_HOME/bin/tcl_setup" \
		>> $(DESTDIR)$(bindir)/orient_mri
	echo "source \$$FREESURFER_HOME/bin/vtk_setup" \
		>> $(DESTDIR)$(bindir)/orient_mri
	echo "source \$$FREESURFER_HOME/bin/kww_setup" \
		>> $(DESTDIR)$(bindir)/orient_mri
	echo "if (-e /bin/pwd) setenv HOME \`cd \$$HOME; /bin/pwd\`" \
		>> $(DESTDIR)$(bindir)/orient_mri
	echo "orient_mri.bin \$$argv" >> $(DESTDIR)$(bindir)/orient_mri
	chmod a+x $(DESTDIR)$(bindir)/orient_mri
	if [[ -f /usr/pubsw/bin/chrpath ]]; then \
	/usr/pubsw/bin/chrpath -d $(DESTDIR)$(bindir)/orient_mri.bin; fi

uninstall-hook:
	rm -f $(DESTDIR)$(bindir)/orient_mri.bin

else
bin_PROGRAMS = 

all:
	@echo "Building KWWidgets apps not enabled."
	@echo "Run configure with either --enable-kww-apps or --with-kwwidgets-dir=DIR"

install-exec-hook:

# these are created on-th-fly (see above), but make dist complains if they
# don't exist, so create some empty files for packaging purposes:
dist-exec-hook:
	@touch vtkKWOrientMRIAppTcl.cxx
	@touch vtkKWOrientMRIWindowTcl.cxx
	@touch vtkKWOrientMRIView2DTcl.cxx
	@touch vtkKWOrientMRIView3DTcl.cxx
	@touch OrientMRILibInit.cxx
	@sleep 2
	@touch vtkKWOrientMRIApp.h
	@touch vtkKWOrientMRIWindow.h
	@touch vtkKWOrientMRIView2D.h
	@touch vtkKWOrientMRIView3D.h
	@touch OrientMRILibInit.data

vtkKWOrientMRIAppTcl.cxx:
vtkKWOrientMRIWindowTcl.cxx:
vtkKWOrientMRIView2DTcl.cxx:
vtkKWOrientMRIView3DTcl.cxx:
OrientMRILibInit.cxx:
endif

clean-local:
	rm -f vtkKWOrientMRIAppTcl.cxx vtkKWOrientMRIWindowTcl.cxx \
	vtkKWOrientMRIView2DTcl.cxx vtkKWOrientMRIView3DTcl.cxx OrientMRILibInit.cxx

EXTRA_DIST=$(orient_mri_DATA) \
OrientMRIEvents.h \
OrientMRILibInit.cxx \
OrientMRILibInit.data \
OrientMRIMain.cxx \
vtkKWOrientMRIApp.cxx \
vtkKWOrientMRIApp.h \
vtkKWOrientMRIAppTcl.cxx \
vtkKWOrientMRIView2D.cxx \
vtkKWOrientMRIView2D.h \
vtkKWOrientMRIView2DTcl.cxx \
vtkKWOrientMRIView3D.cxx \
vtkKWOrientMRIView3D.h \
vtkKWOrientMRIView3DTcl.cxx \
vtkKWOrientMRIWindow.cxx \
vtkKWOrientMRIWindow.h \
vtkKWOrientMRIWindowTcl.cxx \
vtkOrientMRIInteractorStyleView2D.cxx \
vtkOrientMRIInteractorStyleView2D.h \
vtkOrientMRIInteractorStyleView3D.cxx \
vtkOrientMRIInteractorStyleView3D.h


# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra

clean-local:
	-rm *Tcl.cxx OrientMRILibInit.cxx


