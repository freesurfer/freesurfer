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
	ScubaFrame.cpp \
	View.cpp \
	ScubaView.cpp \
	ScubaWindowToRASTranslator.cpp \
	Layer.cpp \
	InputState.cpp \
	ViewState.cpp \
	ScubaToolState.cpp \
	ScubaGlobalPreferences.cpp \
	ProgressDisplayManager.cpp

SCUBA1_SRCS= \
	ToglManager.cpp \
	TclProgressDisplayManager.cpp \
	togl.c

SCUBA2_SRCS= \
	QTScubaWidget.cpp \
	QtProgressDisplayManager.cpp \
	QtScubaUtilities.cpp \
	QtVolumeHistogram.cpp
FORMS = mainwindow.ui \
		viewlayersdlog.ui \
		toolsettingsdlog.ui \
		roisettingsdlog.ui
UICDECLS = .ui/mainwindow.h \
		.ui/viewlayersdlog.h \
		.ui/toolsettingsdlog.h \
		.ui/roisettingsdlog.h
UICIMPLS = .ui/mainwindow.cpp \
		.ui/viewlayersdlog.cpp \
		.ui/toolsettingsdlog.cpp \
		.ui/roisettingsdlog.cpp
SRCMOC   = .moc/moc_QtVolumeHistogram.cpp \
		.moc/moc_QTScubaWidget.cpp \
		.moc/moc_mainwindow.cpp \
		.moc/moc_viewlayersdlog.cpp \
		.moc/moc_toolsettingsdlog.cpp \
		.moc/moc_roisettingsdlog.cpp
OBJMOC = .obj/moc_QtVolumeHistogram.o \
		.obj/moc_QTScubaWidget.o \
		.obj/moc_mainwindow.o \
		.obj/moc_viewlayersdlog.o \
		.obj/moc_toolsettingsdlog.o \
		.obj/moc_roisettingsdlog.o

MOC      = $(QTDIR)/bin/moc
UIC      = $(QTDIR)/bin/uic

bin_PROGRAMS=scuba scuba2
scuba_SOURCES=$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) $(GUI_SRCS) $(LAYER_SRCS) $(SCUBA1_SRCS) scuba.cpp 
LDADD= $(addprefix $(top_srcdir)/,$(LIBS_MGH)) $(LIBS_TCL_OPENGL) 
scuba_LDFLAGS=

scuba2_SOURCES=$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) $(GUI_SRCS) $(LAYER_SRCS) $(MOC_SRCS) $(SCUBA2_SRCS) $(FORMS) $(UICDECLS) $(UICIMPLS) $(SRCMOC) $(OBJMOC) .ui/image_collection.cpp main.cpp
scuba2_LDFLAGS=-L/$(QTDIR)/lib -lqt-mt
scuba2_CXXFLAGS=-DQT_SHARED -DQT_NO_DEBUG -DQT_THREAD_SUPPORT -I$(QTDIR)/include -I.ui/ 

if SCUBA_BUNDLE
scuba2_DEPENDENCIES=scuba2.app
endif


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

# All the scuba2 stuff.
.ui/mainwindow.h: mainwindow.ui
	$(UIC) mainwindow.ui -o .ui/mainwindow.h
.ui/mainwindow.cpp: .ui/mainwindow.h .ui/toolsettingsdlog.h .ui/viewlayersdlog.h .ui/roisettingsdlog.h mainwindow.ui
	$(UIC) mainwindow.ui -i mainwindow.h -o .ui/mainwindow.cpp
.obj/mainwindow.o: .ui/mainwindow.cpp mainwindow.ui.h .ui/mainwindow.h
	$(CXX) $(CXXFLAGS) -o .obj/mainwindow.o .ui/mainwindow.cpp

.ui/toolsettingsdlog.h: toolsettingsdlog.ui
	$(UIC) toolsettingsdlog.ui -o .ui/toolsettingsdlog.h
.ui/toolsettingsdlog.cpp: .ui/mainwindow.h .ui/toolsettingsdlog.h .ui/viewlayersdlog.h .ui/roisettingsdlog.h  toolsettingsdlog.ui
	$(UIC) toolsettingsdlog.ui -i toolsettingsdlog.h -o .ui/toolsettingsdlog.cpp
.obj/toolsettingsdlog.o: .ui/toolsettingsdlog.cpp toolsettingsdlog.ui.h .ui/toolsettingsdlog.h
	$(CXX) $(CXXFLAGS) -o .obj/toolsettingsdlog.o .ui/toolsettingsdlog.cpp

.ui/viewlayersdlog.h: viewlayersdlog.ui
	$(UIC) viewlayersdlog.ui -o .ui/viewlayersdlog.h
.ui/viewlayersdlog.cpp: .ui/mainwindow.h .ui/toolsettingsdlog.h .ui/viewlayersdlog.h .ui/roisettingsdlog.h viewlayersdlog.ui
	$(UIC) viewlayersdlog.ui -i viewlayersdlog.h -o .ui/viewlayersdlog.cpp
.obj/viewlayersdlog.o: .ui/viewlayersdlog.cpp viewlayersdlog.ui.h .ui/viewlayersdlog.h
	$(CXX) $(CXXFLAGS) -o .obj/viewlayersdlog.o .ui/viewlayersdlog.cpp

.ui/roisettingsdlog.h: roisettingsdlog.ui
	$(UIC) roisettingsdlog.ui -o .ui/roisettingsdlog.h
.ui/roisettingsdlog.cpp: .ui/mainwindow.h .ui/toolsettingsdlog.h .ui/viewlayersdlog.h .ui/roisettingsdlog.h roisettingsdlog.ui
	$(UIC) roisettingsdlog.ui -i roisettingsdlog.h -o .ui/roisettingsdlog.cpp
.obj/roisettingsdlog.o: .ui/roisettingsdlog.cpp roisettingsdlog.ui.h .ui/roisettingsdlog.h
	$(CXX) $(CXXFLAGS) -o .obj/roisettingsdlog.o .ui/roisettingsdlog.cpp

.obj/moc_QtVolumeHistogram.o: .moc/moc_QtVolumeHistogram.cpp
	$(CXX) $(CXXFLAGS) -o .obj/moc_QtVolumeHistogram.o
.moc/moc_QtVolumeHistogram.cpp: QtVolumeHistogram.h
	$(MOC) QtVolumeHistogram.h -o .moc/moc_QtVolumeHistogram.cpp

.obj/moc_QTScubaWidget.o: .moc/moc_QTScubaWidget.cpp
	$(CXX) $(CXXFLAGS) -o .obj/moc_QTScubaWidget.o
.moc/moc_QTScubaWidget.cpp: QTScubaWidget.h
	$(MOC) QTScubaWidget.h -o .moc/moc_QTScubaWidget.cpp

.obj/moc_mainwindow.o: .moc/moc_mainwindow.cpp
	$(CXX) $(CXXFLAGS) -o .obj/moc_mainwindow.o
.moc/moc_mainwindow.cpp: .ui/mainwindow.h
	$(MOC) .ui/mainwindow.h -o .moc/moc_mainwindow.cpp

.obj/moc_viewlayersdlog.o: .moc/moc_viewlayersdlog.cpp
	$(CXX) $(CXXFLAGS) -o .obj/moc_viewlayersdlog.o
.moc/moc_viewlayersdlog.cpp: .ui/viewlayersdlog.h
	$(MOC) .ui/viewlayersdlog.h -o .moc/moc_viewlayersdlog.cpp

.obj/moc_toolsettingsdlog.o: .moc/moc_toolsettingsdlog.cpp
	$(CXX) $(CXXFLAGS) -o .obj/moc_toolsettingsdlog.o
.moc/moc_toolsettingsdlog.cpp: .ui/toolsettingsdlog.h
	$(MOC) .ui/toolsettingsdlog.h -o .moc/moc_toolsettingsdlog.cpp

.obj/moc_roisettingsdlog.o: .moc/moc_roisettingsdlog.cpp
	$(CXX) $(CXXFLAGS) -o .obj/moc_roisettingsdlog.o
.moc/moc_roisettingsdlog.cpp: .ui/roisettingsdlog.h
	$(MOC) .ui/roisettingsdlog.h -o .moc/moc_roisettingsdlog.cpp

.obj/image_collection.o: .ui/image_collection.cpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o .obj/image_collection.o .ui/image_collection.cpp

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
		images/icon_folder.png \
		images/icon_home.png \
		images/icon_label_to_marked.png \
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
		images/icon_zoom_in.png \
		images/icon_draw_edge_path.png \
		images/view_inplane.png \
		images/view_inplane_x.png \
		images/view_inplane_y.png \
		icon_eye.xpm \
		icon_lock.xpm \
		icon_pencil.xpm
	$(UIC)  -embed scuba images/filenew images/fileopen images/filesave images/print images/undo images/redo images/editcut images/editcopy images/editpaste images/searchfind images/icon_arrow_ccw.png images/icon_arrow_cw.png images/icon_arrow_down.png images/icon_arrow_expand_x.png images/icon_arrow_expand_y.png images/icon_arrow_left.png images/icon_arrow_right.png images/icon_arrow_rot_x_neg.png images/icon_arrow_rot_x_pos.png images/icon_arrow_rot_y_neg.png images/icon_arrow_rot_y_pos.png images/icon_arrow_rot_z_neg.png images/icon_arrow_rot_z_pos.png images/icon_arrow_shrink_x.png images/icon_arrow_shrink_y.png images/icon_arrow_up.png images/icon_aux_volume.png images/icon_brush_3d.png images/icon_brush_circle.png images/icon_brush_square.png images/icon_color_label.png images/icon_color_picker.png images/icon_control_point.png images/icon_cursor_goto.png images/icon_cursor_save.png images/icon_cut_area.png images/icon_cut_clear.png images/icon_cut_closed_line.png images/icon_cut_line.png images/icon_cut_plane.png images/icon_disk.png images/icon_draw_line.png images/icon_draw_line_closed.png images/icon_edit_ctrlpts.png images/icon_edit_custom.png images/icon_edit_label.png images/icon_edit_parc.png images/icon_edit_volume.png images/icon_erase_label.png images/icon_erase_line.png images/icon_fill.png images/icon_fill_label.png images/icon_folder.png images/icon_home.png images/icon_label_to_marked.png images/icon_line_tool.png images/icon_linked_cursors.png images/icon_main_volume.png images/icon_marked_to_label.png images/icon_marker_crosshair.png images/icon_marker_diamond.png images/icon_navigate.png images/icon_orientation_coronal.png images/icon_orientation_horizontal.png images/icon_orientation_sagittal.png images/icon_redraw.png images/icon_rotate_plane.png images/icon_snapshot_load.png images/icon_snapshot_save.png images/icon_stopwatch.png images/icon_surface_main.png images/icon_surface_original.png images/icon_surface_pial.png images/icon_view_31.png images/icon_view_mosaic.png images/icon_view_multiple.png images/icon_view_single.png images/icon_zoom_in.png images/icon_draw_edge_path.png images/view_inplane.png images/view_inplane_x.png images/view_inplane_y.png icon_eye.xpm icon_lock.xpm icon_pencil.xpm -o .ui/image_collection.cpp


scuba2.app: scuba2.app/Contents/MacOS/../PkgInfo scuba2.app/Contents/MacOS/../Info.plist scuba2.app/Contents/MacOS/scuba2

scuba2.app/Contents/MacOS/../PkgInfo: 
	@test -d scuba2.app/Contents/MacOS/ || mkdir -p scuba2.app/Contents/MacOS/
	@rm -f scuba2.app/Contents/MacOS/../PkgInfo
	@echo "APPL????" >scuba2.app/Contents/MacOS/../PkgInfo

scuba2.app/Contents/MacOS/../Info.plist: 
	@test -d scuba2.app/Contents/MacOS/ || mkdir -p scuba2.app/Contents/MacOS/
	@rm -f scuba2.app/Contents/MacOS/../Info.plist
	@sed -e "s,@ICON@,application.icns,g" -e "s,@EXECUTABLE@,scuba,g" "/space/birn/50/freesurfer/build/storm/qt/mkspecs/default/Info.plist.app" >"scuba2.app/Contents/MacOS/../Info.plist"

scuba2.app/Contents/MacOS/scuba2: scuba2
	@test -d scuba2.app/Contents/MacOS/ || mkdir -p scuba2.app/Contents/MacOS/
	@rm -f scuba2.app/Contents/MacOS/scuba2
	@cp scuba2 scuba2.app/Contents/MacOS/
