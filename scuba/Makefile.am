#
# Makefile.am for scuba
#

AM_CPPFLAGS=-I/usr/include/GL $(GLUT_CFLAGS) $(TCL_CFLAGS) $(tcl_includes) -I$(top_srcdir)/include -I$(top_srcdir)/dicom -I/opt/intel/mkl70/include

AM_LDFLAGS=-L$(top_srcdir)/lib $(LDFLAGS_TCL)

UTILITIES_SRCS=Volume3.cpp \
	Point2.cpp \
	Point3.cpp \
	Path.cpp \
	Array2.cpp \
	Utilities.cpp \
	Matrix44.cpp \
	Transform44.cpp \
	VectorOps.cpp

BASIC_SRCS=TclCommandManager.cpp \
	DebugReporter.cpp \
	IDTracker.cpp \
	PreferencesManager.cpp \
	Broadcaster.cpp \
	Listener.cpp \
	UndoManager.cpp \
	Timer.cpp \
	ShortestPathFinder.cpp \
	PathManager.cpp 

FSDATA_SRCS= \
	DataManager.cpp \
	DataCollection.cpp \
	ScubaDataCollectionFactory.cpp \
	VolumeCollection.cpp \
	SurfaceCollection.cpp \
	ScubaTransform.cpp \
	ScubaColorLUT.cpp \
	ScubaROI.cpp \
	ScubaROIVolume.cpp

LAYER_SRCS=ScubaLayerFactory.cpp \
	ScubaLayer2DMRI.cpp \
	ScubaLayer2DMRIS.cpp

GUI_SRCS= \
	WindowFrame.cpp \
	ToglManager.cpp \
	ScubaFrame.cpp \
	View.cpp \
	ScubaView.cpp \
	ScubaWindowToRASTranslator.cpp \
	Layer.cpp \
	InputState.cpp \
	ViewState.cpp \
	ScubaToolState.cpp \
	ScubaGlobalPreferences.cpp \
	ProgressDisplayManager.cpp \
	TclProgressDisplayManager.cpp \
	togl.c

bin_PROGRAMS = scuba
scuba_SOURCES=$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) $(GUI_SRCS) $(LAYER_SRCS) scuba.cpp 
LDADD= $(addprefix $(top_srcdir)/,$(LIBS_MGH)) $(LIBS_TCL_OPENGL)
scuba_LDFLAGS=

foodir=$(prefix)/lib/tcl
foo_SCRIPTS=scuba.tcl tkcon.tcl histolabel.tcl

check_PROGRAMS= testIDTracker testBroadcaster testUtilities \
	testMatrix44 testTransform44 testVolume3 testPath \
	testShortestPathFinder testToglManager testDataManager \
	testTclCommandManager testPreferencesManager testUndoManager \
	testDataCollection testScubaDataCollections testVolumeCollection \
	testSurfaceCollection testScubaLayers testScubaLayer2DMRI \
	testScubaView testScubaFrame \
	testScubaTransform testScubaColorLUT testScubaROI testScubaROIVolume \
	testVectorOps testRASToVoxel testCoordTimer

#if WANT_QTSCUBA
# bin_PROGRAMS += qtscuba
# qtscuba_SOURCES=$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) $(GUI_SRCS) $(LAYER_SRCS) QTScubaWidget.cpp qtscuba.cpp
# qtscuba_CXXFLAGS=-I${QTDIR}/include -DQT_THREAD_SUPPORT -DQT_SHARED
# qtscuba_LDFLAGS=-L${QTDIR}/lib -lqt-mt
# check_PROGRAMS += testVolumeHistogram
#endif

TESTS=$(check_PROGRAMS) 

## Test programs.
testIDTracker_SOURCES=test_IDTracker.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testBroadcaster_SOURCES=test_Broadcaster.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testUtilities_SOURCES=test_Utilities.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testMatrix44_SOURCES=test_Matrix44.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testTransform44_SOURCES=test_Transform44.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testVolume3_SOURCES=test_Volume3.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testPath_SOURCES=test_Path.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testShortestPathFinder_SOURCES=test_ShortestPathFinder.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testVectorOps_SOURCES=test_VectorOps.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testRASToVoxel_SOURCES= test_RASToVoxel.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS)

testToglManager_SOURCES=test_ToglManager.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(GUI_SRCS) $(FSDATA_SRCS)
testDataManager_SOURCES=test_DataManager.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS)
testTclCommandManager_SOURCES=test_TclCommandManager.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) 
testPreferencesManager_SOURCES=test_PreferencesManager.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) 
testUndoManager_SOURCES=test_UndoManager.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testCoordTimer_SOURCES= test_CoordTimer.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) $(GUI_SRCS)
testCoordTimer_LDFLAGS=-L/opt/intel/mkl70/lib/32 -lmkl -lm

testDataCollection_SOURCES=test_DataCollection.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS)
testScubaDataCollections_SOURCES=test_ScubaDataCollections.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) 
testVolumeCollection_SOURCES=test_VolumeCollection.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) 
testSurfaceCollection_SOURCES= test_SurfaceCollection.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) 

testScubaLayers_SOURCES=test_ScubaLayers.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) $(GUI_SRCS) $(LAYER_SRCS)
testScubaLayer2DMRI_SOURCES=test_ScubaLayer2DMRI.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) $(GUI_SRCS) $(LAYER_SRCS) 

testScubaView_SOURCES=test_ScubaView.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(GUI_SRCS) $(FSDATA_SRCS)
testScubaFrame_SOURCES=test_ScubaFrame.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) $(GUI_SRCS) $(LAYER_SRCS) 
testScubaTransform_SOURCES=test_ScubaTransform.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) 
testScubaColorLUT_SOURCES=test_ScubaColorLUT.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) 
testScubaROI_SOURCES=test_ScubaROI.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) 
testScubaROIVolume_SOURCES=test_ScubaROIVolume.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) 


# if WANT_QTSCUBA
# testVolumeHistogram_SOURCES=$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) $(GUI_SRCS) $(LAYER_SRCS) QtVolumeHistogram.cpp test_VolumeHistogram.cpp
# testVolumeHistogram_CXXFLAGS=-I${QTDIR}/include -DQT_THREAD_SUPPORT -DQT_SHARED
# testVolumeHistogram_LDFLAGS=-L${QTDIR}/lib -lqt-mt
#endif
