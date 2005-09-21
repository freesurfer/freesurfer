#
# Makefile.am for scuba
#

AM_CPPFLAGS=$(GLUT_CFLAGS) \
	$(GL_FLAGS) \
	$(X_CFLAGS) \
	$(TCL_CFLAGS) \
	$(tcl_includes) \
	-I$(top_srcdir)/include \
	-I$(top_srcdir)/dicom \
	-I/opt/intel/mkl70/include

AM_LDFLAGS=-L$(top_srcdir)/lib \
	$(LDFLAGS_TCL) $(GLUT_LIBS)

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
	ScubaROIVolume.cpp \
	SegmentationVolumeReport.cpp

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
	ProgressDisplayManager.cpp \
	ScubaKeyCombo.cpp \
	TclScubaKeyCombo.cpp

SCUBA1_SRCS= \
	ToglManager.cpp \
	TclProgressDisplayManager.cpp \
	togl.c

SCUBA2_SRCS= \
	QTScubaWidget.cpp \
	QtProgressDisplayManager.cpp \
	QtScubaUtilities.cpp \
	QtVolumeHistogram.cpp \
	QtScubaLabelValueArea.cpp \
	QtScubaLineEditWithTclCallback.cpp \
	QtTclListener.cpp

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
		.moc/moc_QtScubaLabelValueArea.cpp \
		.moc/moc_QtScubaLineEditWithTclCallback.cpp \
		.moc/moc_mainwindow.cpp \
		.moc/moc_viewlayersdlog.cpp \
		.moc/moc_toolsettingsdlog.cpp \
		.moc/moc_roisettingsdlog.cpp
OBJMOC = .obj/moc_QtVolumeHistogram.$(OBJEXT) \
		.obj/moc_QTScubaWidget.$(OBJEXT) \
		.obj/moc_QtScubaLabelValueArea.$(OBJEXT) \
		.obj/moc_QtScubaLineEditWithTclCallback.$(OBJEXT) \
		.obj/moc_mainwindow.$(OBJEXT) \
		.obj/moc_viewlayersdlog.$(OBJEXT) \
		.obj/moc_toolsettingsdlog.$(OBJEXT) \
		.obj/moc_roisettingsdlog.$(OBJEXT)

MOC      = $(QTDIR)/bin/moc
UIC      = $(QTDIR)/bin/uic

bin_PROGRAMS=scuba scuba2 ${SCUBA_BUNDLE_TARGET}
EXTRA_PROGRAMS=
scuba_SOURCES=$(UTILITIES_SRCS) \
	$(BASIC_SRCS) \
	$(FSDATA_SRCS) \
	$(GUI_SRCS) \
	$(LAYER_SRCS) \
	$(SCUBA1_SRCS) \
	scuba.cpp 
scuba_LDADD= $(addprefix $(top_srcdir)/,$(LIBS_MGH)) \
	$(LIBS_TCL_OPENGL) $(X_LIBS) -lX11 -lXmu -lGL
scuba_LDFLAGS=

scuba2_SOURCES=$(UTILITIES_SRCS) \
	$(BASIC_SRCS) \
	$(FSDATA_SRCS) \
	$(GUI_SRCS) \
	$(LAYER_SRCS) \
	$(MOC_SRCS) \
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
	testVectorOps testRASToVoxel testCoordTimer \
	testSegmentationVolumeReport testScubaKeyCombo testTclScubaKeyCombo

TESTS=$(check_PROGRAMS) 

## Test programs.
testIDTracker_SOURCES=test_IDTracker.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testIDTracker_LDADD=

testBroadcaster_SOURCES=test_Broadcaster.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testBroadcaster_LDADD=

testUtilities_SOURCES=test_Utilities.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testUtilities_LDADD=

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
testVolume3_LDADD=

testPath_SOURCES=test_Path.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testPath_LDADD=

testShortestPathFinder_SOURCES=test_ShortestPathFinder.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testShortestPathFinder_LDADD=

testVectorOps_SOURCES=test_VectorOps.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testVectorOps_LDADD=

testRASToVoxel_SOURCES= test_RASToVoxel.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS)
testRASToVoxel_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_srcdir)/,$(LIBS_MGH))

testToglManager_SOURCES=test_ToglManager.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(GUI_SRCS) $(FSDATA_SRCS)
testToglManager_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_srcdir)/,$(LIBS_MGH))

testDataManager_SOURCES=test_DataManager.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS)
testDataManager_LDADD=

testTclCommandManager_SOURCES=test_TclCommandManager.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) 
testTclCommandManager_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_srcdir)/,$(LIBS_MGH))

testPreferencesManager_SOURCES=test_PreferencesManager.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) 
testPreferencesManager_LDADD=

testUndoManager_SOURCES=test_UndoManager.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS)
testUndoManager_LDADD=

testCoordTimer_SOURCES= test_CoordTimer.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) $(GUI_SRCS)
testCoordTimer_LDADD=
testCoordTimer_LDFLAGS=-L/opt/intel/mkl70/lib/32 -lmkl -lm

testDataCollection_SOURCES=test_DataCollection.cpp \
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS)
testDataCollection_LDADD=

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

testTclScubaKeyCombo_SOURCES=test_TclScubaKeyCombo.cpp \
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
	$(addprefix $(top_srcdir)/,$(LIBS_MGH))

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
	$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) 
testSegmentationVolumeReport_LDADD=$(LIBS_TCL_OPENGL) \
	$(addprefix $(top_srcdir)/,$(LIBS_MGH))

# All the scuba2 stuff.
.ui/mainwindow.h: ui_dir mainwindow.ui
	$(UIC) mainwindow.ui -o .ui/mainwindow.h
.ui/mainwindow.cpp: ui_dir \
	.ui/mainwindow.h \
	.ui/toolsettingsdlog.h \
	.ui/viewlayersdlog.h \
	.ui/roisettingsdlog.h \
	mainwindow.ui
	$(UIC) mainwindow.ui -i mainwindow.h -o .ui/mainwindow.cpp
.obj/mainwindow.$(OBJEXT): obj_dir .ui/mainwindow.cpp mainwindow.ui.h .ui/mainwindow.h
	$(CXX) $(CXXFLAGS) -o .obj/mainwindow.$(OBJEXT) .ui/mainwindow.cpp

.ui/toolsettingsdlog.h: ui_dir toolsettingsdlog.ui
	$(UIC) toolsettingsdlog.ui -o .ui/toolsettingsdlog.h
.ui/toolsettingsdlog.cpp: ui_dir \
	.ui/mainwindow.h \
	.ui/toolsettingsdlog.h \
	.ui/viewlayersdlog.h \
	.ui/roisettingsdlog.h  \
	toolsettingsdlog.ui
	$(UIC) toolsettingsdlog.ui -i toolsettingsdlog.h -o .ui/toolsettingsdlog.cpp
.obj/toolsettingsdlog.$(OBJEXT): obj_dir \
	.ui/toolsettingsdlog.cpp \
	toolsettingsdlog.ui.h \
	.ui/toolsettingsdlog.h
	$(CXX) $(CXXFLAGS) -o .obj/toolsettingsdlog.$(OBJEXT) .ui/toolsettingsdlog.cpp

.ui/viewlayersdlog.h: ui_dir viewlayersdlog.ui
	$(UIC) viewlayersdlog.ui -o .ui/viewlayersdlog.h
.ui/viewlayersdlog.cpp: ui_dir \
	.ui/mainwindow.h \
	.ui/toolsettingsdlog.h \
	.ui/viewlayersdlog.h \
	.ui/roisettingsdlog.h \
	viewlayersdlog.ui
	$(UIC) viewlayersdlog.ui -i viewlayersdlog.h -o .ui/viewlayersdlog.cpp
.obj/viewlayersdlog.$(OBJEXT): obj_dir \
	.ui/viewlayersdlog.cpp \
	viewlayersdlog.ui.h \
	.ui/viewlayersdlog.h
	$(CXX) $(CXXFLAGS) -o .obj/viewlayersdlog.$(OBJEXT) .ui/viewlayersdlog.cpp

.ui/roisettingsdlog.h: ui_dir roisettingsdlog.ui
	$(UIC) roisettingsdlog.ui -o .ui/roisettingsdlog.h
.ui/roisettingsdlog.cpp: ui_dir \
	.ui/mainwindow.h \
	.ui/toolsettingsdlog.h \
	.ui/viewlayersdlog.h \
	.ui/roisettingsdlog.h \
	roisettingsdlog.ui
	$(UIC) roisettingsdlog.ui -i roisettingsdlog.h -o .ui/roisettingsdlog.cpp
.obj/roisettingsdlog.$(OBJEXT): obj_dir \
	.ui/roisettingsdlog.cpp \
	roisettingsdlog.ui.h \
	.ui/roisettingsdlog.h
	$(CXX) $(CXXFLAGS) -o .obj/roisettingsdlog.$(OBJEXT) .ui/roisettingsdlog.cpp

.obj/moc_QtVolumeHistogram.$(OBJEXT): obj_dir .moc/moc_QtVolumeHistogram.cpp
	$(CXX) $(CXXFLAGS) -o .obj/moc_QtVolumeHistogram.$(OBJEXT)
.moc/moc_QtVolumeHistogram.cpp: moc_dir QtVolumeHistogram.h
	$(MOC) QtVolumeHistogram.h -o .moc/moc_QtVolumeHistogram.cpp

.obj/moc_QTScubaWidget.$(OBJEXT): obj_dir .moc/moc_QTScubaWidget.cpp
	$(CXX) $(CXXFLAGS) -o .obj/moc_QTScubaWidget.$(OBJEXT)
.moc/moc_QTScubaWidget.cpp: moc_dir QTScubaWidget.h
	$(MOC) QTScubaWidget.h -o .moc/moc_QTScubaWidget.cpp

.obj/moc_QtScubaLabelValueArea.$(OBJEXT): obj_dir .moc/moc_QtScubaLabelValueArea.cpp
	$(CXX) $(CXXFLAGS) -o .obj/moc_QtScubaLabelValueArea.$(OBJEXT)
.moc/moc_QtScubaLabelValueArea.cpp: moc_dir QtScubaLabelValueArea.h
	$(MOC) QtScubaLabelValueArea.h -o .moc/moc_QtScubaLabelValueArea.cpp

.obj/moc_QtScubaLineEditWithTclCallback.$(OBJEXT): obj_dir .moc/moc_QtScubaLineEditWithTclCallback.cpp
	$(CXX) $(CXXFLAGS) -o .obj/moc_QtScubaLineEditWithTclCallback.$(OBJEXT)
.moc/moc_QtScubaLineEditWithTclCallback.cpp: moc_dir QtScubaLineEditWithTclCallback.h
	$(MOC) QtScubaLineEditWithTclCallback.h -o .moc/moc_QtScubaLineEditWithTclCallback.cpp

.obj/moc_mainwindow.$(OBJEXT): obj_dir .moc/moc_mainwindow.cpp
	$(CXX) $(CXXFLAGS) -o .obj/moc_mainwindow.$(OBJEXT)
.moc/moc_mainwindow.cpp: moc_dir .ui/mainwindow.h
	$(MOC) .ui/mainwindow.h -o .moc/moc_mainwindow.cpp

.obj/moc_viewlayersdlog.$(OBJEXT): obj_dir .moc/moc_viewlayersdlog.cpp
	$(CXX) $(CXXFLAGS) -o .obj/moc_viewlayersdlog.$(OBJEXT)
.moc/moc_viewlayersdlog.cpp: moc_dir .ui/viewlayersdlog.h
	$(MOC) .ui/viewlayersdlog.h -o .moc/moc_viewlayersdlog.cpp

.obj/moc_toolsettingsdlog.$(OBJEXT): obj_dir .moc/moc_toolsettingsdlog.cpp
	$(CXX) $(CXXFLAGS) -o .obj/moc_toolsettingsdlog.$(OBJEXT)
.moc/moc_toolsettingsdlog.cpp: moc_dir .ui/toolsettingsdlog.h
	$(MOC) .ui/toolsettingsdlog.h -o .moc/moc_toolsettingsdlog.cpp

.obj/moc_roisettingsdlog.$(OBJEXT): obj_dir .moc/moc_roisettingsdlog.cpp
	$(CXX) $(CXXFLAGS) -o .obj/moc_roisettingsdlog.$(OBJEXT)
.moc/moc_roisettingsdlog.cpp: moc_dir .ui/roisettingsdlog.h
	$(MOC) .ui/roisettingsdlog.h -o .moc/moc_roisettingsdlog.cpp

.obj/image_collection.$(OBJEXT): obj_dir .ui/image_collection.cpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o .obj/image_collection.$(OBJEXT) \
	.ui/image_collection.cpp

.ui/image_collection.cpp: ui_dir images/filenew \
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
	$(UIC)  -embed scuba2 images/filenew images/fileopen images/filesave images/print images/undo images/redo images/editcut images/editcopy images/editpaste images/searchfind images/icon_arrow_ccw.png images/icon_arrow_cw.png images/icon_arrow_down.png images/icon_arrow_expand_x.png images/icon_arrow_expand_y.png images/icon_arrow_left.png images/icon_arrow_right.png images/icon_arrow_rot_x_neg.png images/icon_arrow_rot_x_pos.png images/icon_arrow_rot_y_neg.png images/icon_arrow_rot_y_pos.png images/icon_arrow_rot_z_neg.png images/icon_arrow_rot_z_pos.png images/icon_arrow_shrink_x.png images/icon_arrow_shrink_y.png images/icon_arrow_up.png images/icon_aux_volume.png images/icon_brush_3d.png images/icon_brush_circle.png images/icon_brush_square.png images/icon_color_label.png images/icon_color_picker.png images/icon_control_point.png images/icon_cursor_goto.png images/icon_cursor_save.png images/icon_cut_area.png images/icon_cut_clear.png images/icon_cut_closed_line.png images/icon_cut_line.png images/icon_cut_plane.png images/icon_disk.png images/icon_draw_line.png images/icon_draw_line_closed.png images/icon_edit_ctrlpts.png images/icon_edit_custom.png images/icon_edit_label.png images/icon_edit_parc.png images/icon_edit_volume.png images/icon_erase_label.png images/icon_erase_line.png images/icon_fill.png images/icon_fill_label.png images/icon_folder.png images/icon_home.png images/icon_label_to_marked.png images/icon_line_tool.png images/icon_linked_cursors.png images/icon_main_volume.png images/icon_marked_to_label.png images/icon_marker_crosshair.png images/icon_marker_diamond.png images/icon_navigate.png images/icon_orientation_coronal.png images/icon_orientation_horizontal.png images/icon_orientation_sagittal.png images/icon_redraw.png images/icon_rotate_plane.png images/icon_snapshot_load.png images/icon_snapshot_save.png images/icon_stopwatch.png images/icon_surface_main.png images/icon_surface_original.png images/icon_surface_pial.png images/icon_view_31.png images/icon_view_mosaic.png images/icon_view_multiple.png images/icon_view_single.png images/icon_zoom_in.png images/icon_draw_edge_path.png images/view_inplane.png images/view_inplane_x.png images/view_inplane_y.png icon_eye.xpm icon_lock.xpm icon_pencil.xpm images/icon_zoom_out.png images/icon_fill_roi.png images/icon_fill_volume.png images/icon_label_list.png images/icon_label_table.png -o .ui/image_collection.cpp


scuba_bundle: scuba2
	@test -d scuba2.app/Contents/MacOS/ || mkdir -p scuba2.app/Contents/MacOS/
	@rm -f scuba2.app/Contents/MacOS/../PkgInfo
	@echo "APPL????" >scuba2.app/Contents/MacOS/../PkgInfo
	@test -d scuba2.app/Contents/MacOS/ || mkdir -p scuba2.app/Contents/MacOS/
	@rm -f scuba2.app/Contents/MacOS/../Info.plist
	@sed -e "s,@ICON@,application.icns,g" -e "s,@EXECUTABLE@,scuba,g" "$(QTDIR)/mkspecs/default/Info.plist.app" >"scuba2.app/Contents/MacOS/../Info.plist"
	@test -d scuba2.app/Contents/MacOS/ || mkdir -p scuba2.app/Contents/MacOS/
	@rm -f scuba2.app/Contents/MacOS/scuba2
	@cp scuba2 scuba2.app/Contents/MacOS/
	@cp -r scuba2.app $(DESTDIR)$(bindir)

obj_dir:
	@test -d .obj || mkdir -p .obj

moc_dir:
	@test -d .moc || mkdir -p .moc

ui_dir:
	@test -d .ui || mkdir -p .ui

.PHONY: scuba_bundle ui_dir obj_dir moc_dir

