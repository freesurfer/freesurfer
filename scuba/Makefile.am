#
# Makefile.am for scuba
#

if ENABLE_QT_APPS
	SUBDIRS=images
endif

AM_CPPFLAGS=$(GLUT_CFLAGS) \
	$(GL_FLAGS) \
	$(X_CFLAGS) \
	$(TCL_CFLAGS) \
	$(tcl_includes) \
	-I$(top_srcdir)/include \
	-I$(top_srcdir)/dicom \
	-I/opt/intel/mkl70/include \
	-Iqconsole -Iqconsole/include

AM_LDFLAGS=$(LDFLAGS_TCL) $(GLUT_LIBS)

UTILITIES_SRCS=Volume3.cpp Volume3.h \
	Point2.cpp Point2.h \
	Point3.cpp Point3.h \
	Path.cpp Path.h \
	Array2.cpp Array2.h \
	Utilities.cpp Utilities.h \
	Matrix44.cpp Matrix44.h \
	Transform44.cpp Transform44.h \
	VectorOps.cpp VectorOps.h \
	string_fixed.h togl.h

BASIC_SRCS=TclCommandManager.cpp TclCommandManager.h \
	DebugReporter.cpp DebugReporter.h \
	IDTracker.cpp IDTracker.h \
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

QCONSOLE_SRCS = \
	qconsole/src/commands.cpp qconsole/include/commands.h \
	qconsole/src/commandsManager.cpp qconsole/include/commandsManager.h \
	qconsole/src/interceptor.cpp qconsole/include/interceptor.h \
	qconsole/src/qconsole.cpp qconsole/include/qconsole.h \
	qconsole/src/qtclconsole.cpp qconsole/include/qtclconsole.h

SCUBA1_SRCS= \
	ToglManager.cpp ToglManager.h \
	TclProgressDisplayManager.cpp TclProgressDisplayManager.h \
	TclChartWindow.cpp TclChartWindow.h \
	togl.c \
	TclScubaKeyCombo.cpp TclScubaKeyCombo.h

if ENABLE_QT_APPS
SCUBA2_SRCS= \
	QTScubaWidget.cpp QTScubaWidget.h \
	QtProgressDisplayManager.cpp QtProgressDisplayManager.h \
	QtScubaUtilities.cpp QtScubaUtilities.h \
	QtVolumeHistogram.cpp QtVolumeHistogram.h \
	QtScubaLabelValueArea.cpp QtScubaLabelValueArea.h \
	QtScubaLineEditWithTclCallback.cpp QtScubaLineEditWithTclCallback.h \
	QtTclListener.cpp QtTclListener.h \
	QtScubaKeyCombo.cpp QtScubaKeyCombo.h \
	QtKeyGrabberButton.cpp QtKeyGrabberButton.h

FORMS = mainwindow.ui \
		viewlayersdlog.ui \
		toolsettingsdlog.ui \
		roisettingsdlog.ui \
		prefsdlog.ui \
		editdatdlog.ui \
		findsurfacevertexdlog.ui
UICDECLS = .ui/mainwindow.h mainwindow.ui.h \
		.ui/viewlayersdlog.h viewlayersdlog.ui.h \
		.ui/toolsettingsdlog.h toolsettingsdlog.ui.h \
		.ui/roisettingsdlog.h roisettingsdlog.ui.h \
		.ui/prefsdlog.h prefsdlog.ui.h \
		.ui/editdatdlog.h editdatdlog.ui.h \
		.ui/findsurfacevertexdlog.h findsurfacevertexdlog.ui.h

UICIMPLS = .ui/mainwindow.cpp \
		.ui/viewlayersdlog.cpp \
		.ui/toolsettingsdlog.cpp \
		.ui/roisettingsdlog.cpp \
		.ui/prefsdlog.cpp \
		.ui/editdatdlog.cpp \
		.ui/findsurfacevertexdlog.cpp
SRCMOC   = .moc/moc_QtVolumeHistogram.cpp \
		.moc/moc_QTScubaWidget.cpp \
		.moc/moc_QtScubaLabelValueArea.cpp \
		.moc/moc_QtScubaLineEditWithTclCallback.cpp \
		.moc/moc_QtKeyGrabberButton.cpp \
		.moc/moc_mainwindow.cpp \
		.moc/moc_viewlayersdlog.cpp \
		.moc/moc_toolsettingsdlog.cpp \
		.moc/moc_roisettingsdlog.cpp \
		.moc/moc_prefsdlog.cpp \
		.moc/moc_editdatdlog.cpp \
		.moc/moc_findsurfacevertexdlog.cpp \
		.moc/moc_interceptor.cpp \
		.moc/moc_qconsole.cpp

#OBJMOC = moc_QtVolumeHistogram.$(OBJEXT) \
#		moc_QTScubaWidget.$(OBJEXT) \
#		moc_QtScubaLabelValueArea.$(OBJEXT) \
#		moc_QtScubaLineEditWithTclCallback.$(OBJEXT) \
#		moc_mainwindow.$(OBJEXT) \
#		moc_viewlayersdlog.$(OBJEXT) \
#		moc_toolsettingsdlog.$(OBJEXT) \
#		moc_roisettingsdlog.$(OBJEXT)

MOC      = $(QTDIR)/bin/moc
UIC      = $(QTDIR)/bin/uic

# endif ENABLE_QT_APPS:
endif

if ENABLE_TCLTK_APPS
  SCUBA_TCLTK=scuba
endif
if ENABLE_QT_APPS
  SCUBA_QT=scuba2 ${SCUBA_BUNDLE_TARGET}
endif
bin_PROGRAMS=${SCUBA_TCLTK} ${SCUBA_QT}

EXTRA_PROGRAMS=
scuba_SOURCES=$(UTILITIES_SRCS) \
	$(BASIC_SRCS) \
	$(FSDATA_SRCS) \
	$(GUI_SRCS) \
	$(LAYER_SRCS) \
	$(CHART_SRCS) \
	$(SCUBA1_SRCS) \
	scuba.cpp
scuba_LDADD= $(addprefix $(top_srcdir)/,$(LIBS_MGH)) \
	$(LIBS_TCL_OPENGL) $(X_LIBS) -lX11 -lXmu -lGL
scuba_LDFLAGS=

if ENABLE_QT_APPS
scuba2_SOURCES=$(UTILITIES_SRCS) \
	$(BASIC_SRCS) \
	$(FSDATA_SRCS) \
	$(CHART_SRCS) \
	$(GUI_SRCS) \
	$(LAYER_SRCS) \
	$(MOC_SRCS) \
	$(QCONSOLE_SRCS) \
	$(SCUBA2_SRCS) \
	$(FORMS) \
	$(UICDECLS) \
	$(UICIMPLS) \
	$(SRCMOC) \
	$(OBJMOC) \
	.ui/image_collection.cpp \
	main.cpp
scuba2_LDADD= $(addprefix $(top_srcdir)/,$(LIBS_MGH)) \
	$(LIBS_TCL_OPENGL) \
	-L$(QTDIR)/lib -lqt-mt
scuba2_LDFLAGS=
scuba2_CXXFLAGS=-DQT_SHARED \
	-DQT_NO_DEBUG \
	-DQT_THREAD_SUPPORT \
	-I$(QTDIR)/include \
	-I.ui/
endif

foodir=$(prefix)/lib/tcl
foo_SCRIPTS=scuba.tcl tkcon.tcl histolabel.tcl TclChartWindow.tcl

check_PROGRAMS= testIDTracker testBroadcaster testUtilities \
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

EXTRA_DIST=$(wildcard *.xpm) $(foo_SCRIPTS) $(wildcard *.tcl) \
	testSetup.csh testCleanup.csh test_data.tar.gz \
	testSegmentationVolumeReport-correct.txt \
	testSegmentationVolumeReportIntReport-correct.txt \
	testLUT.lut

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
	$(addprefix $(top_srcdir)/,$(LIBS_MGH))

testBroadcaster_SOURCES=test_Broadcaster.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testBroadcaster_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_srcdir)/,$(LIBS_MGH))

testUtilities_SOURCES=test_Utilities.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testUtilities_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_srcdir)/,$(LIBS_MGH))

testMatrix44_SOURCES=test_Matrix44.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testMatrix44_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_srcdir)/,$(LIBS_MGH))

testTransform44_SOURCES=test_Transform44.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testTransform44_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_srcdir)/,$(LIBS_MGH))

testVolume3_SOURCES=test_Volume3.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testVolume3_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_srcdir)/,$(LIBS_MGH))

testPath_SOURCES=test_Path.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testPath_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_srcdir)/,$(LIBS_MGH))

testShortestPathFinder_SOURCES=test_ShortestPathFinder.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testShortestPathFinder_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_srcdir)/,$(LIBS_MGH))

testVectorOps_SOURCES=test_VectorOps.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testVectorOps_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_srcdir)/,$(LIBS_MGH))

testRASToVoxel_SOURCES= test_RASToVoxel.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS)
testRASToVoxel_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_srcdir)/,$(LIBS_MGH))

testToglManager_SOURCES=test_ToglManager.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(GUI_SRCS) $(FSDATA_SRCS)
testToglManager_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_srcdir)/,$(LIBS_MGH)) $(GLUT_LIBS)

testDataManager_SOURCES=test_DataManager.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS)
testDataManager_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_srcdir)/,$(LIBS_MGH))

testTclCommandManager_SOURCES=test_TclCommandManager.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testTclCommandManager_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_srcdir)/,$(LIBS_MGH))

testPreferencesManager_SOURCES=test_PreferencesManager.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testPreferencesManager_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_srcdir)/,$(LIBS_MGH))

testUndoManager_SOURCES=test_UndoManager.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testUndoManager_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_srcdir)/,$(LIBS_MGH))

testCoordTimer_SOURCES= test_CoordTimer.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) $(GUI_SRCS)
testCoordTimer_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_srcdir)/,$(LIBS_MGH))

testDataCollection_SOURCES=test_DataCollection.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS)
testDataCollection_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_srcdir)/,$(LIBS_MGH))

testScubaDataCollections_SOURCES=test_ScubaDataCollections.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS)
testScubaDataCollections_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_srcdir)/,$(LIBS_MGH))

testVolumeCollection_SOURCES=test_VolumeCollection.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS)
testVolumeCollection_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_srcdir)/,$(LIBS_MGH))

testSurfaceCollection_SOURCES= test_SurfaceCollection.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS)
testSurfaceCollection_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_srcdir)/,$(LIBS_MGH))

testScubaKeyCombo_SOURCES=test_ScubaKeyCombo.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) $(GUI_SRCS)
testScubaKeyCombo_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_srcdir)/,$(LIBS_MGH))

testTclScubaKeyCombo_SOURCES=test_TclScubaKeyCombo.cpp TclScubaKeyCombo.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) $(GUI_SRCS)
testTclScubaKeyCombo_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_srcdir)/,$(LIBS_MGH))

testScubaLayers_SOURCES=test_ScubaLayers.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) $(GUI_SRCS) $(LAYER_SRCS)
testScubaLayers_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_srcdir)/,$(LIBS_MGH))

testScubaLayer2DMRI_SOURCES=test_ScubaLayer2DMRI.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) $(GUI_SRCS) $(LAYER_SRCS)
testScubaLayer2DMRI_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_srcdir)/,$(LIBS_MGH))

testScubaView_SOURCES=test_ScubaView.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(GUI_SRCS) $(FSDATA_SRCS) $(SCUBA1_SRCS)
testScubaView_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_srcdir)/,$(LIBS_MGH)) $(X_LIBS)

testScubaFrame_SOURCES=test_ScubaFrame.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) $(GUI_SRCS) $(LAYER_SRCS)
testScubaFrame_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_srcdir)/,$(LIBS_MGH))

testScubaTransform_SOURCES=test_ScubaTransform.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS)
testScubaTransform_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_srcdir)/,$(LIBS_MGH))

testScubaColorLUT_SOURCES=test_ScubaColorLUT.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS)
testScubaColorLUT_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_srcdir)/,$(LIBS_MGH))

testScubaROI_SOURCES=test_ScubaROI.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS)
testScubaROI_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_srcdir)/,$(LIBS_MGH))

testScubaROIVolume_SOURCES=test_ScubaROIVolume.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS)
testScubaROIVolume_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_srcdir)/,$(LIBS_MGH))

testSegmentationVolumeReport_SOURCES=test_SegmentationVolumeReport.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) $(GUI_SRCS) $(LAYER_SRCS) $(CHART_SRCS)
testSegmentationVolumeReport_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_srcdir)/,$(LIBS_MGH))

testTclChartWindow_SOURCES=test_TclChartWindow.cpp TclChartWindow.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) $(GUI_SRCS) $(LAYER_SRCS) $(CHART_SRCS)
testTclChartWindow_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_srcdir)/,$(LIBS_MGH))

if ENABLE_QT_APPS
# All the scuba2 stuff.
.ui/mainwindow.h: mainwindow.ui
	$(UIC) mainwindow.ui -o $@
.ui/mainwindow.cpp: \
	.ui/mainwindow.h \
	.ui/toolsettingsdlog.h \
	.ui/viewlayersdlog.h \
	.ui/roisettingsdlog.h \
	.ui/prefsdlog.h \
	.ui/editdatdlog.h \
	.ui/findsurfacevertexdlog.h \
	mainwindow.ui
	$(UIC) mainwindow.ui -i mainwindow.h -o $@
mainwindow.$(OBJEXT): .ui/mainwindow.cpp
	mainwindow.ui.h .ui/mainwindow.h

.ui/toolsettingsdlog.h: toolsettingsdlog.ui
	$(UIC) toolsettingsdlog.ui -o $@
.ui/toolsettingsdlog.cpp: \
	.ui/mainwindow.h \
	.ui/toolsettingsdlog.h \
	.ui/viewlayersdlog.h \
	.ui/roisettingsdlog.h  \
	.ui/prefsdlog.h \
	.ui/editdatdlog.h \
	.ui/findsurfacevertexdlog.h \
	toolsettingsdlog.ui
	$(UIC) toolsettingsdlog.ui -i toolsettingsdlog.h -o $@
toolsettingsdlog.$(OBJEXT): .ui/toolsettingsdlog.cpp
	toolsettingsdlog.ui.h .ui/toolsettingsdlog.h

.ui/viewlayersdlog.h: viewlayersdlog.ui
	$(UIC) viewlayersdlog.ui -o $@
.ui/viewlayersdlog.cpp: \
	.ui/mainwindow.h \
	.ui/toolsettingsdlog.h \
	.ui/viewlayersdlog.h \
	.ui/roisettingsdlog.h \
	.ui/prefsdlog.h \
	.ui/editdatdlog.h \
	.ui/findsurfacevertexdlog.h \
	viewlayersdlog.ui
	$(UIC) viewlayersdlog.ui -i viewlayersdlog.h -o $@
viewlayersdlog.$(OBJEXT): .ui/viewlayersdlog.cpp \
	viewlayersdlog.ui.h .ui/viewlayersdlog.h

.ui/roisettingsdlog.h: roisettingsdlog.ui
	$(UIC) roisettingsdlog.ui -o $@
.ui/roisettingsdlog.cpp: \
	.ui/mainwindow.h \
	.ui/toolsettingsdlog.h \
	.ui/viewlayersdlog.h \
	.ui/roisettingsdlog.h \
	.ui/prefsdlog.h \
	.ui/editdatdlog.h \
	.ui/findsurfacevertexdlog.h \
	roisettingsdlog.ui
	$(UIC) roisettingsdlog.ui -i roisettingsdlog.h -o $@
roisettingsdlog.$(OBJEXT): .ui/roisettingsdlog.cpp \
	roisettingsdlog.ui.h .ui/roisettingsdlog.h

.ui/prefsdlog.h: prefsdlog.ui
	$(UIC) prefsdlog.ui -o $@
.ui/prefsdlog.cpp: \
	.ui/mainwindow.h \
	.ui/toolsettingsdlog.h \
	.ui/viewlayersdlog.h \
	.ui/roisettingsdlog.h \
	.ui/prefsdlog.h \
	.ui/editdatdlog.h \
	.ui/findsurfacevertexdlog.h \
	prefsdlog.ui
	$(UIC) prefsdlog.ui -i prefsdlog.h -o $@
prefsdlog.$(OBJEXT): .ui/prefsdlog.cpp \
	prefsdlog.ui.h .ui/prefsdlog.h

.ui/editdatdlog.h: editdatdlog.ui
	$(UIC) editdatdlog.ui -o $@
.ui/editdatdlog.cpp: \
	.ui/mainwindow.h \
	.ui/toolsettingsdlog.h \
	.ui/viewlayersdlog.h \
	.ui/roisettingsdlog.h \
	.ui/prefsdlog.h \
	.ui/editdatdlog.h \
	.ui/findsurfacevertexdlog.h \
	editdatdlog.ui
	$(UIC) editdatdlog.ui -i editdatdlog.h -o $@
editdatdlog.$(OBJEXT): .ui/editdatdlog.cpp \
	editdatdlog.ui.h .ui/editdatdlog.h

.ui/findsurfacevertexdlog.h: findsurfacevertexdlog.ui
	$(UIC) findsurfacevertexdlog.ui -o $@
.ui/findsurfacevertexdlog.cpp: \
	.ui/mainwindow.h \
	.ui/toolsettingsdlog.h \
	.ui/viewlayersdlog.h \
	.ui/roisettingsdlog.h \
	.ui/prefsdlog.h \
	.ui/editdatdlog.h \
	.ui/findsurfacevertexdlog.h \
	findsurfacevertexdlog.ui
	$(UIC) findsurfacevertexdlog.ui -i findsurfacevertexdlog.h -o $@
findsurfacevertexdlog.$(OBJEXT): .ui/findsurfacevertexdlog.cpp \
	findsurfacevertexdlog.ui.h .ui/findsurfacevertexdlog.h

moc_QtVolumeHistogram.$(OBJEXT): .moc/moc_QtVolumeHistogram.cpp
.moc/moc_QtVolumeHistogram.cpp: QtVolumeHistogram.h
	$(MOC) QtVolumeHistogram.h -o $@

moc_QTScubaWidget.$(OBJEXT): .moc/moc_QTScubaWidget.cpp
.moc/moc_QTScubaWidget.cpp: QTScubaWidget.h
	$(MOC) QTScubaWidget.h -o $@

moc_QtScubaLabelValueArea.$(OBJEXT): .moc/moc_QtScubaLabelValueArea.cpp
.moc/moc_QtScubaLabelValueArea.cpp: QtScubaLabelValueArea.h
	$(MOC) QtScubaLabelValueArea.h -o $@

moc_QtScubaLineEditWithTclCallback.$(OBJEXT): .moc/moc_QtScubaLineEditWithTclCallback.cpp
.moc/moc_QtScubaLineEditWithTclCallback.cpp: QtScubaLineEditWithTclCallback.h
	$(MOC) QtScubaLineEditWithTclCallback.h -o $@

moc_QtKeyGrabberButton.$(OBJEXT): .moc/moc_QtKeyGrabberButton.cpp
.moc/moc_QtKeyGrabberButton.cpp: QtKeyGrabberButton.h
	$(MOC) QtKeyGrabberButton.h -o $@

moc_mainwindow.$(OBJEXT): .moc/moc_mainwindow.cpp
.moc/moc_mainwindow.cpp: .ui/mainwindow.h
	$(MOC) .ui/mainwindow.h -o $@

moc_viewlayersdlog.$(OBJEXT): .moc/moc_viewlayersdlog.cpp
.moc/moc_viewlayersdlog.cpp: .ui/viewlayersdlog.h
	$(MOC) .ui/viewlayersdlog.h -o $@

moc_toolsettingsdlog.$(OBJEXT): .moc/moc_toolsettingsdlog.cpp
.moc/moc_toolsettingsdlog.cpp: .ui/toolsettingsdlog.h
	$(MOC) .ui/toolsettingsdlog.h -o $@

moc_roisettingsdlog.$(OBJEXT): .moc/moc_roisettingsdlog.cpp
.moc/moc_roisettingsdlog.cpp: .ui/roisettingsdlog.h
	$(MOC) .ui/roisettingsdlog.h -o $@

moc_prefsdlog.$(OBJEXT): .moc/moc_prefsdlog.cpp
.moc/moc_prefsdlog.cpp: .ui/prefsdlog.h
	$(MOC) .ui/prefsdlog.h -o $@

moc_editdatdlog.$(OBJEXT): .moc/moc_editdatdlog.cpp
.moc/moc_editdatdlog.cpp: .ui/editdatdlog.h
	$(MOC) .ui/editdatdlog.h -o $@

moc_findsurfacevertexdlog.$(OBJEXT): .moc/moc_findsurfacevertexdlog.cpp
.moc/moc_findsurfacevertexdlog.cpp: .ui/findsurfacevertexdlog.h
	$(MOC) .ui/findsurfacevertexdlog.h -o $@

moc_interceptor.$(OBJEXT): .moc/moc_interceptor.cpp
.moc/moc_interceptor.cpp:
	$(MOC) qconsole/include/interceptor.h -o $@

moc_qconsole.$(OBJEXT): .moc/moc_qconsole.cpp
.moc/moc_qconsole.cpp:
	$(MOC) qconsole/include/qconsole.h -o $@

image_collection.$(OBJEXT): .ui/image_collection.cpp

.ui/image_collection.cpp: images/filenew \
		images/fileopen \
		images/filesave \
		images/print \
		images/undo \
		images/redo \
		images/editcut \
		images/editcopy \
		images/editpaste \
		images/searchfind \
		images/icon_arrow_ccw.png \
		images/icon_arrow_cw.png \
		images/icon_arrow_down.png \
		images/icon_arrow_expand_x.png \
		images/icon_arrow_expand_y.png \
		images/icon_arrow_left.png \
		images/icon_arrow_right.png \
		images/icon_arrow_rot_x_neg.png \
		images/icon_arrow_rot_x_pos.png \
		images/icon_arrow_rot_y_neg.png \
		images/icon_arrow_rot_y_pos.png \
		images/icon_arrow_rot_z_neg.png \
		images/icon_arrow_rot_z_pos.png \
		images/icon_arrow_shrink_x.png \
		images/icon_arrow_shrink_y.png \
		images/icon_arrow_up.png \
		images/icon_aux_volume.png \
		images/icon_brush_3d.png \
		images/icon_brush_circle.png \
		images/icon_brush_square.png \
		images/icon_color_label.png \
		images/icon_color_picker.png \
		images/icon_control_point.png \
		images/icon_cursor_goto.png \
		images/icon_cursor_save.png \
		images/icon_cut_area.png \
		images/icon_cut_clear.png \
		images/icon_cut_closed_line.png \
		images/icon_cut_line.png \
		images/icon_cut_plane.png \
		images/icon_disk.png \
		images/icon_draw_line.png \
		images/icon_draw_line_closed.png \
		images/icon_edit_ctrlpts.png \
		images/icon_edit_custom.png \
		images/icon_edit_label.png \
		images/icon_edit_parc.png \
		images/icon_edit_volume.png \
		images/icon_erase_label.png \
		images/icon_erase_line.png \
		images/icon_fill.png \
		images/icon_fill_label.png \
		images/icon_fill_roi.png \
		images/icon_fill_volume.png \
		images/icon_folder.png \
		images/icon_home.png \
		images/icon_label_to_marked.png \
		images/icon_label_table.png \
		images/icon_label_list.png \
		images/icon_line_tool.png \
		images/icon_linked_cursors.png \
		images/icon_main_volume.png \
		images/icon_marked_to_label.png \
		images/icon_marker_crosshair.png \
		images/icon_marker_diamond.png \
		images/icon_navigate.png \
		images/icon_orientation_coronal.png \
		images/icon_orientation_horizontal.png \
		images/icon_orientation_sagittal.png \
		images/icon_redraw.png \
		images/icon_rotate_plane.png \
		images/icon_snapshot_load.png \
		images/icon_snapshot_save.png \
		images/icon_stopwatch.png \
		images/icon_surface_main.png \
		images/icon_surface_original.png \
		images/icon_surface_pial.png \
		images/icon_view_31.png \
		images/icon_view_mosaic.png \
		images/icon_view_multiple.png \
		images/icon_view_single.png \
		images/icon_zoom_out.png \
		images/icon_zoom_in.png \
		images/icon_draw_edge_path.png \
		images/view_inplane.png \
		images/view_inplane_x.png \
		images/view_inplane_y.png \
		icon_eye.xpm \
		icon_lock.xpm \
		icon_pencil.xpm
	$(UIC)  -embed scuba2 images/filenew images/fileopen images/filesave images/print images/undo images/redo images/editcut images/editcopy images/editpaste images/searchfind images/icon_arrow_ccw.png images/icon_arrow_cw.png images/icon_arrow_down.png images/icon_arrow_expand_x.png images/icon_arrow_expand_y.png images/icon_arrow_left.png images/icon_arrow_right.png images/icon_arrow_rot_x_neg.png images/icon_arrow_rot_x_pos.png images/icon_arrow_rot_y_neg.png images/icon_arrow_rot_y_pos.png images/icon_arrow_rot_z_neg.png images/icon_arrow_rot_z_pos.png images/icon_arrow_shrink_x.png images/icon_arrow_shrink_y.png images/icon_arrow_up.png images/icon_aux_volume.png images/icon_brush_3d.png images/icon_brush_circle.png images/icon_brush_square.png images/icon_color_label.png images/icon_color_picker.png images/icon_control_point.png images/icon_cursor_goto.png images/icon_cursor_save.png images/icon_cut_area.png images/icon_cut_clear.png images/icon_cut_closed_line.png images/icon_cut_line.png images/icon_cut_plane.png images/icon_disk.png images/icon_draw_line.png images/icon_draw_line_closed.png images/icon_edit_ctrlpts.png images/icon_edit_custom.png images/icon_edit_label.png images/icon_edit_parc.png images/icon_edit_volume.png images/icon_erase_label.png images/icon_erase_line.png images/icon_fill.png images/icon_fill_label.png images/icon_folder.png images/icon_home.png images/icon_label_to_marked.png images/icon_line_tool.png images/icon_linked_cursors.png images/icon_main_volume.png images/icon_marked_to_label.png images/icon_marker_crosshair.png images/icon_marker_diamond.png images/icon_navigate.png images/icon_orientation_coronal.png images/icon_orientation_horizontal.png images/icon_orientation_sagittal.png images/icon_redraw.png images/icon_rotate_plane.png images/icon_snapshot_load.png images/icon_snapshot_save.png images/icon_stopwatch.png images/icon_surface_main.png images/icon_surface_original.png images/icon_surface_pial.png images/icon_view_31.png images/icon_view_mosaic.png images/icon_view_multiple.png images/icon_view_single.png images/icon_zoom_in.png images/icon_draw_edge_path.png images/view_inplane.png images/view_inplane_x.png images/view_inplane_y.png icon_eye.xpm icon_lock.xpm icon_pencil.xpm images/icon_zoom_out.png images/icon_fill_roi.png images/icon_fill_volume.png images/icon_label_list.png images/icon_label_table.png -o $@

scuba_bundle: scuba2
	@echo Making OSX bundle...
	@test -d scuba2.app/Contents/MacOS/ || mkdir -p scuba2.app/Contents/MacOS/
	@rm -f scuba2.app/Contents/MacOS/../PkgInfo
	@echo "APPL????" >scuba2.app/Contents/MacOS/../PkgInfo
	@test -d scuba2.app/Contents/MacOS/ || mkdir -p scuba2.app/Contents/MacOS/
	@rm -f scuba2.app/Contents/MacOS/../Info.plist
	@sed -e "s,@ICON@,application.icns,g" -e "s,@EXECUTABLE@,scuba2,g" "$(QTDIR)/mkspecs/default/Info.plist.app" >"scuba2.app/Contents/MacOS/../Info.plist"
	@test -d scuba2.app/Contents/MacOS/ || mkdir -p scuba2.app/Contents/MacOS/
	@rm -f scuba2.app/Contents/MacOS/scuba2
	@cp scuba2 scuba2.app/Contents/MacOS/
	@test -d scuba2.app/Contents/Frameworks || mkdir -p scuba2.app/Contents/Frameworks
	@cp $(QTDIR)/lib/libqt-mt.3.dylib scuba2.app/Contents/Frameworks
	@install_name_tool -id @executable_path/../Frameworks/libqt-mt.3.dylib scuba2.app/Contents/Frameworks/libqt-mt.3.dylib
	@install_name_tool -change libqt-mt.3.dylib @executable_path/../Frameworks/libqt-mt.3.dylib scuba2.app/Contents/MacOS/scuba2
	@cp -r scuba2.app $(DESTDIR)$(bindir)
	@echo "$(FREESURFER_HOME)/bin/scuba2.app/Contents/MacOS/scuba2" > $(DESTDIR)$(bindir)/scuba2
	@chmod a+x $(DESTDIR)$(bindir)/scuba2

# endif ENABLE_QT_APPS:
endif

.PHONY: scuba_bundle

# put a wrapper around the bin, used to setup tcltktixblt enviro vars
install-exec-hook:
	mv $(bindir)/scuba $(bindir)/scuba.bin
	echo "#!/bin/tcsh -ef" > $(bindir)/scuba
	echo "source \$$FREESURFER_HOME/bin/tcl_setup" >> $(bindir)/scuba
	echo "scuba.bin \$$argv" >> $(bindir)/scuba
	chmod a+x $(bindir)/scuba

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
