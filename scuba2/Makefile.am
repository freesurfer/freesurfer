if BUILD_SCUBA

AM_CPPFLAGS=-I$(top_srcdir)/include

AM_CXXFLAGS=-Wno-deprecated -I$(top_srcdir)/include \
	-I$(top_srcdir)/vtkfsio \
	-I$(top_srcdir)/vtkutils \
	$(VTK_CXXFLAGS) $(KWWIDGETS_CXXFLAGS) $(TCL_CFLAGS)

if ENABLE_KWWIDGETS_APPS
bin_PROGRAMS = scuba2
scuba2_SOURCES= \
	ScubaMain.cxx string_fixed.h \
	ScubaViewProperties.h \
	vtkKWScubaApp.cxx vtkKWScubaApp.h vtkKWScubaAppTcl.cxx \
	vtkKWScubaApplicationSettingsInterface.cxx \
	vtkKWScubaApplicationSettingsInterface.h \
	vtkKWScubaApplicationSettingsInterfaceTcl.cxx \
	vtkKWScubaWindow.cxx vtkKWScubaWindow.h vtkKWScubaWindowTcl.cxx \
	vtkKWScubaView.cxx vtkKWScubaView.h vtkKWScubaViewTcl.cxx \
	vtkScubaInteractorStyle.cxx vtkScubaInteractorStyle.h \
	vtkScubaInteractorStyleTcl.cxx \
	vtkKWScubaTool.cxx vtkKWScubaTool.h vtkKWScubaToolTcl.cxx \
	vtkKWScubaToolNavigate.cxx vtkKWScubaToolNavigate.h \
	vtkKWScubaToolNavigateTcl.cxx \
	vtkKWScubaToolEdit2DMRI.cxx vtkKWScubaToolEdit2DMRI.h \
	vtkKWScubaToolEdit2DMRITcl.cxx \
	ScubaInfoItem.cxx ScubaInfoItem.h \
	vtkKWScubaLayer.cxx vtkKWScubaLayer.h vtkKWScubaLayerTcl.cxx\
	vtkKWScubaLayerCollection.cxx vtkKWScubaLayerCollection.h \
	vtkKWScubaLayerCollectionTcl.cxx\
	vtkKWScubaLayerCollectionMRI.cxx vtkKWScubaLayerCollectionMRI.h \
	vtkKWScubaLayerCollectionMRITcl.cxx\
	vtkKWScubaLayerCollectionMRIS.cxx vtkKWScubaLayerCollectionMRIS.h \
	vtkKWScubaLayerCollectionMRISTcl.cxx\
	vtkKWScubaLayerCollectionDTI.cxx vtkKWScubaLayerCollectionDTI.h \
	vtkKWScubaLayerCollectionDTITcl.cxx\
	vtkKWScubaLayerCollectionPath.cxx vtkKWScubaLayerCollectionPath.h \
	vtkKWScubaLayerCollectionPathTcl.cxx\
	vtkKWScubaLayerCollectionODF.cxx vtkKWScubaLayerCollectionODF.h \
	vtkKWScubaLayerCollectionODFTcl.cxx\
	ScubaCollectionProperties.h \
	ScubaCollectionPropertiesMRI.h \
	ScubaCollectionPropertiesMRIS.h \
	ScubaCollectionPropertiesDTI.h \
	ScubaCollectionPropertiesODF.h \
	ScubaCollectionPropertiesPath.h \
	vtkKWScubaLayer2DMRI.cxx vtkKWScubaLayer2DMRI.h \
	vtkKWScubaLayer2DMRITcl.cxx \
	vtkKWScubaLayer3DMRI.cxx vtkKWScubaLayer3DMRI.h \
	vtkKWScubaLayer3DMRITcl.cxx \
	vtkKWScubaLayer2DMRIS.cxx vtkKWScubaLayer2DMRIS.h \
	vtkKWScubaLayer2DMRISTcl.cxx \
	vtkKWScubaLayer3DMRIS.cxx vtkKWScubaLayer3DMRIS.h \
	vtkKWScubaLayer3DMRISTcl.cxx \
	vtkKWScubaLayer2DDTI.cxx vtkKWScubaLayer2DDTI.h \
	vtkKWScubaLayer2DDTITcl.cxx \
	vtkKWScubaLayer3DDTI.cxx vtkKWScubaLayer3DDTI.h \
	vtkKWScubaLayer3DDTITcl.cxx \
	vtkKWScubaLayer2DPath.cxx vtkKWScubaLayer2DPath.h \
	vtkKWScubaLayer2DPathTcl.cxx \
	vtkKWScubaLayer3DPath.cxx vtkKWScubaLayer3DPath.h \
	vtkKWScubaLayer3DPathTcl.cxx \
	vtkKWScubaLayer3DODF.cxx vtkKWScubaLayer3DODF.h \
	vtkKWScubaLayer3DODFTcl.cxx \
	Broadcaster.cxx Broadcaster.h Listener.cxx Listener.h \
	IDTracker.cxx IDTracker.h \
	ScubaLibInit.cxx

scuba2_LDADD=  \
	$(top_builddir)/fsgdf/libfsgdf.a \
	$(top_builddir)/vtkfsio/libvtkfsio.a \
	$(top_builddir)/vtkutils/libvtkutils.a \
	$(addprefix $(top_builddir)/, $(LIBS_MGH)) $(LIBS_TCL_OPENGL) \
	$(VTK_LIBS) $(VTK_VERDICT_LIB) -lvtkWidgets -lvtkGraphics \
	-lvtkImaging -lvtkFiltering \
	-lvtkCommon -lvtksys \
	-lvtkHybrid -lvtkexoIIc -lvtkNetCDF \
	-lvtkVolumeRendering -lvtkRendering -lvtkftgl \
	$(VTK_EXPAT_LIB) $(VTK_FREETYPE_LIB) \
	$(VTK_TIFF_LIB) $(VTK_JPEG_LIB) $(VTK_SQLITE_LIB) \
	$(VTK_PNG_LIB) $(VTK_Z_LIB) $(VTK_METAIO_LIB)\
	-lvtkIO	-lvtkDICOMParser $(VTK_MPEG2ENCODE_LIB) \
	-lvtkImagingTCL -lvtkVolumeRenderingTCL -lvtkRenderingTCL \
	-lvtkFilteringTCL -lvtkWidgetsTCL \
	-lvtkHybridTCL -lvtkGraphicsTCL -lvtkImagingTCL -lvtkIOTCL \
	-lvtkCommonTCL \
	$(KWWIDGETS_LIBS) -lKWWidgets
scuba2_LDFLAGS=$(OS_LDFLAGS) $(LDFLAGS_TCL)

scuba2dir=$(prefix)/lib/resource
scuba2_DATA=ScubaIcons.txt

vtkKWScubaAppTcl.cxx: vtkKWScubaApp.h
	${VTK_DIR}/bin/vtkWrapTcl vtkKWScubaApp.h \
	${VTK_LIB_VTK}/hints 1 vtkKWScubaAppTcl.cxx

vtkKWScubaApplicationSettingsInterfaceTcl.cxx: vtkKWScubaApplicationSettingsInterface.h
	${VTK_DIR}/bin/vtkWrapTcl vtkKWScubaApplicationSettingsInterface.h \
	${VTK_LIB_VTK}/hints 1 vtkKWScubaApplicationSettingsInterfaceTcl.cxx

vtkKWScubaWindowTcl.cxx: vtkKWScubaWindow.h
	${VTK_DIR}/bin/vtkWrapTcl vtkKWScubaWindow.h \
	${VTK_LIB_VTK}/hints 1 vtkKWScubaWindowTcl.cxx

vtkKWScubaViewTcl.cxx: vtkKWScubaView.h
	${VTK_DIR}/bin/vtkWrapTcl vtkKWScubaView.h \
	${VTK_LIB_VTK}/hints 1 vtkKWScubaViewTcl.cxx

vtkKWScubaLayerTcl.cxx: vtkKWScubaLayer.h
	${VTK_DIR}/bin/vtkWrapTcl vtkKWScubaLayer.h \
	${VTK_LIB_VTK}/hints 1 vtkKWScubaLayerTcl.cxx

vtkKWScubaLayerCollectionTcl.cxx: vtkKWScubaLayerCollection.h
	${VTK_DIR}/bin/vtkWrapTcl vtkKWScubaLayerCollection.h \
	${VTK_LIB_VTK}/hints 1 vtkKWScubaLayerCollectionTcl.cxx

vtkKWScubaLayerCollectionMRITcl.cxx: vtkKWScubaLayerCollectionMRI.h
	${VTK_DIR}/bin/vtkWrapTcl vtkKWScubaLayerCollectionMRI.h \
	${VTK_LIB_VTK}/hints 1 vtkKWScubaLayerCollectionMRITcl.cxx

vtkKWScubaLayerCollectionMRISTcl.cxx: vtkKWScubaLayerCollectionMRIS.h
	${VTK_DIR}/bin/vtkWrapTcl vtkKWScubaLayerCollectionMRIS.h \
	${VTK_LIB_VTK}/hints 1 vtkKWScubaLayerCollectionMRISTcl.cxx

vtkKWScubaLayerCollectionDTITcl.cxx: vtkKWScubaLayerCollectionDTI.h
	${VTK_DIR}/bin/vtkWrapTcl vtkKWScubaLayerCollectionDTI.h \
	${VTK_LIB_VTK}/hints 1 vtkKWScubaLayerCollectionDTITcl.cxx

vtkKWScubaLayerCollectionPathTcl.cxx: vtkKWScubaLayerCollectionPath.h
	${VTK_DIR}/bin/vtkWrapTcl vtkKWScubaLayerCollectionPath.h \
	${VTK_LIB_VTK}/hints 1 vtkKWScubaLayerCollectionPathTcl.cxx

vtkKWScubaLayerCollectionODFTcl.cxx: vtkKWScubaLayerCollectionODF.h
	${VTK_DIR}/bin/vtkWrapTcl vtkKWScubaLayerCollectionODF.h \
	${VTK_LIB_VTK}/hints 1 vtkKWScubaLayerCollectionODFTcl.cxx

vtkKWScubaLayer2DMRITcl.cxx: vtkKWScubaLayer2DMRI.h
	${VTK_DIR}/bin/vtkWrapTcl vtkKWScubaLayer2DMRI.h \
	${VTK_LIB_VTK}/hints 1 vtkKWScubaLayer2DMRITcl.cxx

vtkKWScubaLayer3DMRITcl.cxx: vtkKWScubaLayer3DMRI.h
	${VTK_DIR}/bin/vtkWrapTcl vtkKWScubaLayer3DMRI.h \
	${VTK_LIB_VTK}/hints 1 vtkKWScubaLayer3DMRITcl.cxx

vtkKWScubaLayer2DMRISTcl.cxx: vtkKWScubaLayer2DMRIS.h
	${VTK_DIR}/bin/vtkWrapTcl vtkKWScubaLayer2DMRIS.h \
	${VTK_LIB_VTK}/hints 1 vtkKWScubaLayer2DMRISTcl.cxx

vtkKWScubaLayer3DMRISTcl.cxx: vtkKWScubaLayer3DMRIS.h
	${VTK_DIR}/bin/vtkWrapTcl vtkKWScubaLayer3DMRIS.h \
	${VTK_LIB_VTK}/hints 1 vtkKWScubaLayer3DMRISTcl.cxx

vtkKWScubaLayer2DDTITcl.cxx: vtkKWScubaLayer2DDTI.h
	${VTK_DIR}/bin/vtkWrapTcl vtkKWScubaLayer2DDTI.h \
	${VTK_LIB_VTK}/hints 1 vtkKWScubaLayer2DDTITcl.cxx

vtkKWScubaLayer3DDTITcl.cxx: vtkKWScubaLayer3DDTI.h
	${VTK_DIR}/bin/vtkWrapTcl vtkKWScubaLayer3DDTI.h \
	${VTK_LIB_VTK}/hints 1 vtkKWScubaLayer3DDTITcl.cxx

vtkKWScubaLayer2DPathTcl.cxx: vtkKWScubaLayer2DPath.h
	${VTK_DIR}/bin/vtkWrapTcl vtkKWScubaLayer2DPath.h \
	${VTK_LIB_VTK}/hints 1 vtkKWScubaLayer2DPathTcl.cxx

vtkKWScubaLayer3DPathTcl.cxx: vtkKWScubaLayer3DPath.h
	${VTK_DIR}/bin/vtkWrapTcl vtkKWScubaLayer3DPath.h \
	${VTK_LIB_VTK}/hints 1 vtkKWScubaLayer3DPathTcl.cxx

vtkKWScubaLayer3DODFTcl.cxx: vtkKWScubaLayer3DODF.h
	${VTK_DIR}/bin/vtkWrapTcl vtkKWScubaLayer3DODF.h \
	${VTK_LIB_VTK}/hints 1 vtkKWScubaLayer3DODFTcl.cxx

vtkScubaInteractorStyleTcl.cxx: vtkScubaInteractorStyle.h
	${VTK_DIR}/bin/vtkWrapTcl vtkScubaInteractorStyle.h \
	${VTK_LIB_VTK}/hints 1 vtkScubaInteractorStyleTcl.cxx

vtkKWScubaToolTcl.cxx: vtkKWScubaTool.h
	${VTK_DIR}/bin/vtkWrapTcl vtkKWScubaTool.h \
	${VTK_LIB_VTK}/hints 1 vtkKWScubaToolTcl.cxx

vtkKWScubaToolNavigateTcl.cxx: vtkKWScubaToolNavigate.h
	${VTK_DIR}/bin/vtkWrapTcl vtkKWScubaToolNavigate.h \
	${VTK_LIB_VTK}/hints 1 vtkKWScubaToolNavigateTcl.cxx

vtkKWScubaToolEdit2DMRITcl.cxx: vtkKWScubaToolEdit2DMRI.h
	${VTK_DIR}/bin/vtkWrapTcl vtkKWScubaToolEdit2DMRI.h \
	${VTK_LIB_VTK}/hints 1 vtkKWScubaToolEdit2DMRITcl.cxx

ScubaLibInit.cxx: ScubaLibInit.data
	${VTK_DIR}/bin/vtkWrapTclInit ScubaLibInit.data ScubaLibInit.cxx



# HACK - njs: as of 8 Aug 2007, this fails compilation with:
#   vtkKWScubaView.cxx:261: internal compiler error: in extract_insn
# source of error unknown, but disappears if compiled with -O1 instead of -O2
vtkKWScubaView.o: vtkKWScubaView.cxx
	@echo "compiling vtkKWScubaView.cxx to vtkKWScubaView.o (with -O1)"
	source='$<' object='$@' libtool=no \
	depfile='$(DEPDIR)/$*.Po' tmpdepfile='$(DEPDIR)/$*.TPo' \
	$(CXXDEPMODE) $(depcomp) \
	$(CXXCOMPILE) -O1 -c `test -f '$<' || echo '$(srcdir)/'`$<


# put a wrapper around the bin, used to setup tcltktixblt,vtk,kww enviro vars
install-exec-hook:scuba2
	cp $(top_builddir)/scuba2/scuba2 $(DESTDIR)$(bindir)/scuba2.bin
	echo "#!/bin/tcsh -ef" > $(DESTDIR)$(bindir)/scuba2
	echo "source \$$FREESURFER_HOME/bin/tcl_setup" \
		>> $(DESTDIR)$(bindir)/scuba2
	echo "source \$$FREESURFER_HOME/bin/vtk_setup" \
		>> $(DESTDIR)$(bindir)/scuba2
	echo "source \$$FREESURFER_HOME/bin/kww_setup" \
		>> $(DESTDIR)$(bindir)/scuba2
	echo "if (-e /bin/pwd) setenv HOME \`cd \$$HOME; /bin/pwd\`" \
		>> $(DESTDIR)$(bindir)/scuba2
	echo "scuba2.bin \$$argv" >> $(DESTDIR)$(bindir)/scuba2
	chmod a+x $(DESTDIR)$(bindir)/scuba2
	if [[ -f /usr/pubsw/bin/chrpath ]]; then \
	/usr/pubsw/bin/chrpath -d $(DESTDIR)$(bindir)/scuba2.bin; fi

uninstall-hook:
	rm -f $(DESTDIR)$(bindir)/scuba2.bin

clean-local:
	-rm *Tcl.cxx

else
bin_PROGRAMS = 

all:
	@echo "Building KWWidgets apps not enabled."
	@echo "Run configure with either --enable-kww-apps or --with-kwwidgets-dir=DIR"

install-exec-hook:

# these are created on-th-fly (see above), but make dist complains if they
# don't exist, so create some empty files for packaging purposes:
dist-exec-hook:
	@touch vtkKWRGBATransferFunctionEditorTcl.cxx
	@touch vtkKWScubaAppTcl.cxx
	@touch vtkKWScubaLayer2DDTITcl.cxx
	@touch vtkKWScubaLayer2DMRISTcl.cxx
	@touch vtkKWScubaLayer2DMRITcl.cxx
	@touch vtkKWScubaLayer2DPathTcl.cxx
	@touch vtkKWScubaLayer3DDTITcl.cxx
	@touch vtkKWScubaLayer3DMRISTcl.cxx
	@touch vtkKWScubaLayer3DMRITcl.cxx
	@touch vtkKWScubaLayer3DPathTcl.cxx
	@touch vtkKWScubaLayer3DODFTcl.cxx
	@touch vtkKWScubaLayerCollectionDTITcl.cxx
	@touch vtkKWScubaLayerCollectionMRISTcl.cxx
	@touch vtkKWScubaLayerCollectionMRITcl.cxx
	@touch vtkKWScubaLayerCollectionPathTcl.cxx
	@touch vtkKWScubaLayerCollectionODFTcl.cxx
	@touch vtkKWScubaLayerCollectionTcl.cxx
	@touch vtkKWScubaLayerTcl.cxx
	@touch vtkKWScubaToolEdit2DMRITcl.cxx
	@touch vtkKWScubaToolNavigateTcl.cxx
	@touch vtkKWScubaToolTcl.cxx
	@touch vtkKWScubaViewTcl.cxx
	@touch vtkKWScubaWindowTcl.cxx
	@touch ScubaLibInit.cxx
	@sleep 2
	@touch ScubaLibInit.data

vtkKWScubaAppTcl.cxx:
vtkKWScubaWindowTcl.cxx:
vtkKWScubaViewTcl.cxx:
vtkKWScubaLayerTcl.cxx:
vtkKWScubaLayerCollectionTcl.cxx:
vtkKWScubaLayerCollectionMRITcl.cxx:
vtkKWScubaLayerCollectionMRISTcl.cxx:
vtkKWScubaLayerCollectionDTITcl.cxx:
vtkKWScubaLayerCollectionPathTcl.cxx:
vtkKWScubaLayerCollectionODFTcl.cxx:
vtkKWScubaLayer2DMRITcl.cxx:
vtkKWScubaLayer3DMRITcl.cxx:
vtkKWScubaLayer2DMRISTcl.cxx:
vtkKWScubaLayer3DMRISTcl.cxx:
vtkKWScubaLayer2DDTITcl.cxx:
vtkKWScubaLayer3DDTITcl.cxx:
vtkKWScubaLayer2DPathTcl.cxx:
vtkKWScubaLayer3DPathTcl.cxx:
vtkKWScubaLayer3DODFTcl.cxx:
vtkScubaInteractorStyleTcl.cxx:
vtkKWScubaToolTcl.cxx:
vtkKWScubaToolNavigateTcl.cxx:
vtkKWScubaToolEdit2DMRITcl.cxx:
ScubaLibInit.cxx:
endif

EXTRA_DIST=ScubaIcons.txt ScubaIcons.txt ScubaLibInit.data \
	$(scuba2_SOURCES)

# BUILD_SCUBA
endif

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra

