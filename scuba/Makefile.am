#
# Makefile.am for scuba
#

if BUILD_SCUBA

if USE_LOCAL_TIFF
TIFFCFLAGS=-I$(top_srcdir)/tiff
LIBTIFF=$(top_builddir)/tiff/libtiff.a
else
TIFFCFLAGS=$(TIFF_CFLAGS)
LIBTIFF=$(LIB_TIFF)
endif

AM_CPPFLAGS=\
	$(GL_CFLAGS) $(TIFFCFLAGS) \
	$(TCL_CFLAGS) \
	$(tcl_includes) \
	-I$(top_srcdir)/include \
	-I$(top_srcdir)/dicom \
	-I$(top_srcdir)/glut

# -Wno-non-virtual-dtor suppresses valid scuba and cppunit warnings:
# "has virtual functions but non-virtual destructor"
# these warnings appear with gcc 4.#
AM_CXXFLAGS=-Wno-non-virtual-dtor

AM_LDFLAGS=$(LDFLAGS_TCL) $(TIFF_LIBS) $(top_builddir)/glut/libglut.a

UTILITIES_SRCS=Volume3.cpp Volume3.h \
	Point2.h \
	Point3.h \
	Path.cpp Path.h \
	Array2.cpp Array2.h \
	Utilities.cpp Utilities.h \
	Matrix44.cpp Matrix44.h \
	Transform44.cpp Transform44.h \
	VectorOps.cpp VectorOps.h \
	string_fixed.h togl.h

BASIC_SRCS=TclCommandManager.cpp TclCommandManager.h \
	DebugReporter.cpp DebugReporter.h \
	IDTracker.h \
	PreferencesManager.cpp PreferencesManager.h \
	Broadcaster.cpp Broadcaster.h \
	Listener.cpp Listener.h \
	UndoManager.cpp UndoManager.h \
	Timer.cpp Timer.h \
	ShortestPathFinder.cpp ShortestPathFinder.h \
	PathManager.cpp PathManager.h \
	Scuba-impl.h

FSDATA_SRCS= \
	DataManager.cpp DataManager.h \
	DataCollection.cpp DataCollection.h \
	ScubaDataCollectionFactory.cpp ScubaDataCollectionFactory.h \
	VolumeCollection.cpp VolumeCollection.h \
	SurfaceCollection.cpp SurfaceCollection.h \
	ScubaTransform.cpp ScubaTransform.h \
	ScubaColorLUT.cpp ScubaColorLUT.h \
	ScubaROI.cpp ScubaROI.h \
	ScubaROIVolume.cpp ScubaROIVolume.h

LAYER_SRCS=ScubaLayerFactory.cpp ScubaLayerFactory.h \
	ScubaLayer2DMRI.cpp ScubaLayer2DMRI.h \
	ScubaLayer2DMRIS.cpp ScubaLayer2DMRIS.h

GUI_SRCS= \
	WindowFrame.cpp WindowFrame.h \
	ScubaFrame.cpp ScubaFrame.h \
	View.cpp View.h \
	ScubaView.cpp ScubaView.h \
	ScubaWindowToRASTranslator.cpp ScubaWindowToRASTranslator.h \
	Layer.cpp Layer.h \
	InputState.cpp InputState.h \
	ViewState.cpp ViewState.h \
	ScubaToolState.cpp ScubaToolState.h \
	ScubaGlobalPreferences.cpp ScubaGlobalPreferences.h \
	ProgressDisplayManager.cpp ProgressDisplayManager.h \
	ScubaKeyCombo.cpp ScubaKeyCombo.h \
	SegmentationVolumeReport.cpp SegmentationVolumeReport.h \
	ChartWindow.cpp ChartWindow.h

CHART_SRCS= \
	ScubaVolumeROIIntensityChart.cpp \
	ScubaVolumeROIIntensityChart.h \
	ScubaMultiFrameVolumeChart.cpp \
	ScubaMultiFrameVolumeChart.h

SCUBA1_SRCS= \
	ToglManager.cpp ToglManager.h \
	TclProgressDisplayManager.cpp TclProgressDisplayManager.h \
	TclChartWindow.cpp TclChartWindow.h \
	togl.c \
	TclScubaKeyCombo.cpp TclScubaKeyCombo.h

if ENABLE_TCLTK_APPS
  SCUBA_TCLTK=scuba
else
all:
	@echo "Building Tcl/Tk apps not enabled."
	@echo "Run configure with either --enable-tcltk-apps or --with-tcl-dir=DIR"
endif
bin_PROGRAMS=${SCUBA_TCLTK}

EXTRA_PROGRAMS=
scuba_SOURCES=$(UTILITIES_SRCS) \
	$(BASIC_SRCS) \
	$(FSDATA_SRCS) \
	$(GUI_SRCS) \
	$(LAYER_SRCS) \
	$(CHART_SRCS) \
	$(SCUBA1_SRCS) \
	scuba.cpp
scuba_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH)) \
	$(LIBTIFF) $(LIBS_TCL_OPENGL) $(X_LIBS) -lX11 -lXmu -lGL
scuba_LDFLAGS=$(AM_LDFLAGS)

foodir=$(prefix)/lib/tcl
foo_SCRIPTS=\
	histolabel.tcl \
	scuba_glut.tcl \
	scuba.tcl \
	TclChartWindow.tcl \
	test_layers.tcl \
	test_ScubaFrame.tcl \
	test_ScubaLayers.tcl \
	test_TclCommandManager.tcl \
	test_tkUtils.tcl \
	test_ToglManager.tcl \
	test_view.tcl \
	tkcon.tcl

if ENABLE_TCLTK_APPS
# trick to get test data into the build directory
foo_SOURCES=
foo:
	if ! test -f "$(top_builddir)/scuba/test_data.tar.gz"; then \
		cp -v $(top_srcdir)/scuba/test_data.tar.gz \
					$(top_srcdir)/scuba/test_TclCommandManager.tcl \
					$(top_srcdir)/scuba/testLUT.lut \
					$(top_srcdir)/scuba/testSegmentationVolumeReport-correct.txt \
				$(top_srcdir)/scuba/testSegmentationVolumeReportIntReport-correct.txt \
					$(top_builddir)/scuba/; fi;

check_PROGRAMS=foo testIDTracker testBroadcaster testUtilities \
	testMatrix44 testTransform44 testVolume3 testPath \
	testShortestPathFinder testToglManager testDataManager \
	testTclCommandManager testPreferencesManager testUndoManager \
	testDataCollection testScubaDataCollections \
	testSurfaceCollection testScubaLayers testScubaLayer2DMRI \
	testScubaView testScubaFrame \
	testScubaTransform testScubaColorLUT testScubaROI testScubaROIVolume \
	testVectorOps testRASToVoxel testCoordTimer \
	testSegmentationVolumeReport testScubaKeyCombo testTclScubaKeyCombo \
	testVolumeCollection testTclChartWindow
endif

EXTRA_DIST=$(foo_SCRIPTS) \
	testSetup.csh testCleanup.csh test_data.tar.gz \
	testSegmentationVolumeReport-correct.txt \
	testSegmentationVolumeReportIntReport-correct.txt \
	testLUT.lut test_TclCommandManager.tcl

if ENABLE_TCLTK_APPS
TESTS=testSetup.csh \
	testIDTracker testBroadcaster testUtilities \
	testMatrix44 testTransform44 testVolume3 testPath \
	testShortestPathFinder testToglManager testDataManager \
	testTclCommandManager testPreferencesManager \
	testUndoManager \
	testDataCollection testScubaDataCollections \
	testVolumeCollection \
	testSurfaceCollection testScubaLayers testScubaLayer2DMRI \
	testScubaFrame \
	testScubaTransform testScubaColorLUT testScubaROI \
	testScubaROIVolume \
	testVectorOps testRASToVoxel testCoordTimer \
	testSegmentationVolumeReport testScubaKeyCombo \
	testTclScubaKeyCombo \
	testCleanup.csh

## Test programs.
testIDTracker_SOURCES=test_IDTracker.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testIDTracker_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_builddir)/, $(LIBS_MGH))

testBroadcaster_SOURCES=test_Broadcaster.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testBroadcaster_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_builddir)/, $(LIBS_MGH))

testUtilities_SOURCES=test_Utilities.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testUtilities_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_builddir)/, $(LIBS_MGH))

testMatrix44_SOURCES=test_Matrix44.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testMatrix44_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_builddir)/, $(LIBS_MGH))

testTransform44_SOURCES=test_Transform44.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testTransform44_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_builddir)/, $(LIBS_MGH))

testVolume3_SOURCES=test_Volume3.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testVolume3_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_builddir)/, $(LIBS_MGH))

testPath_SOURCES=test_Path.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testPath_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_builddir)/, $(LIBS_MGH))

testShortestPathFinder_SOURCES=test_ShortestPathFinder.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testShortestPathFinder_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_builddir)/, $(LIBS_MGH))

testVectorOps_SOURCES=test_VectorOps.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testVectorOps_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_builddir)/, $(LIBS_MGH))

testRASToVoxel_SOURCES= test_RASToVoxel.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS)
testRASToVoxel_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_builddir)/, $(LIBS_MGH))

testToglManager_SOURCES=test_ToglManager.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(GUI_SRCS) $(FSDATA_SRCS)
testToglManager_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_builddir)/, $(LIBS_MGH)) $(GLUT_LIBS)

testDataManager_SOURCES=test_DataManager.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS)
testDataManager_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_builddir)/, $(LIBS_MGH))

testTclCommandManager_SOURCES=test_TclCommandManager.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testTclCommandManager_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_builddir)/, $(LIBS_MGH))

testPreferencesManager_SOURCES=test_PreferencesManager.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testPreferencesManager_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_builddir)/, $(LIBS_MGH))

testUndoManager_SOURCES=test_UndoManager.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testUndoManager_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_builddir)/, $(LIBS_MGH))

testCoordTimer_SOURCES= test_CoordTimer.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) $(GUI_SRCS)
testCoordTimer_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_builddir)/, $(LIBS_MGH))

testDataCollection_SOURCES=test_DataCollection.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS)
testDataCollection_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_builddir)/, $(LIBS_MGH))

testScubaDataCollections_SOURCES=test_ScubaDataCollections.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS)
testScubaDataCollections_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_builddir)/, $(LIBS_MGH))

testVolumeCollection_SOURCES=test_VolumeCollection.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS)
testVolumeCollection_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_builddir)/, $(LIBS_MGH))

testSurfaceCollection_SOURCES= test_SurfaceCollection.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS)
testSurfaceCollection_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_builddir)/, $(LIBS_MGH))

testScubaKeyCombo_SOURCES=test_ScubaKeyCombo.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) $(GUI_SRCS)
testScubaKeyCombo_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_builddir)/, $(LIBS_MGH))

testTclScubaKeyCombo_SOURCES=test_TclScubaKeyCombo.cpp TclScubaKeyCombo.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) $(GUI_SRCS)
testTclScubaKeyCombo_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_builddir)/, $(LIBS_MGH))

testScubaLayers_SOURCES=test_ScubaLayers.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) $(GUI_SRCS) $(LAYER_SRCS)
testScubaLayers_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_builddir)/, $(LIBS_MGH))

testScubaLayer2DMRI_SOURCES=test_ScubaLayer2DMRI.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) $(GUI_SRCS) $(LAYER_SRCS)
testScubaLayer2DMRI_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_builddir)/, $(LIBS_MGH))

testScubaView_SOURCES=test_ScubaView.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(GUI_SRCS) $(FSDATA_SRCS) $(SCUBA1_SRCS)
testScubaView_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_builddir)/, $(LIBS_MGH)) $(X_LIBS)

testScubaFrame_SOURCES=test_ScubaFrame.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) $(GUI_SRCS) $(LAYER_SRCS)
testScubaFrame_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_builddir)/, $(LIBS_MGH))

testScubaTransform_SOURCES=test_ScubaTransform.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS)
testScubaTransform_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_builddir)/, $(LIBS_MGH))

testScubaColorLUT_SOURCES=test_ScubaColorLUT.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS)
testScubaColorLUT_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_builddir)/, $(LIBS_MGH))

testScubaROI_SOURCES=test_ScubaROI.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS)
testScubaROI_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_builddir)/, $(LIBS_MGH))

testScubaROIVolume_SOURCES=test_ScubaROIVolume.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS)
testScubaROIVolume_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_builddir)/, $(LIBS_MGH))

testSegmentationVolumeReport_SOURCES=test_SegmentationVolumeReport.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) $(GUI_SRCS) \
	$(LAYER_SRCS) $(CHART_SRCS)
testSegmentationVolumeReport_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_builddir)/, $(LIBS_MGH))

testTclChartWindow_SOURCES=test_TclChartWindow.cpp TclChartWindow.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) $(GUI_SRCS) \
	$(LAYER_SRCS) $(CHART_SRCS)
testTclChartWindow_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_builddir)/, $(LIBS_MGH))
endif

if ENABLE_TCLTK_APPS
# put a wrapper around the bin, used to setup tcltktixblt enviro vars
install-exec-hook:
	cp $(top_builddir)/scuba/scuba $(DESTDIR)$(bindir)/scuba.bin
	echo "#!/bin/tcsh -ef" > $(DESTDIR)$(bindir)/scuba
	echo "source \$$FREESURFER_HOME/bin/tcl_setup" >> $(DESTDIR)$(bindir)/scuba
	echo "scuba.bin \$$argv" >> $(DESTDIR)$(bindir)/scuba
	chmod a+x $(DESTDIR)$(bindir)/scuba
	if [[ -f /usr/pubsw/bin/chrpath ]]; then \
	/usr/pubsw/bin/chrpath -d $(DESTDIR)$(bindir)/scuba.bin; fi

uninstall-hook:
	rm -f $(DESTDIR)$(bindir)/scuba.bin

else
install-exec-hook:
endif

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra

# BUILD_SCUBA:
endif
