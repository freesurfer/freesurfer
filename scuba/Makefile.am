#
# Makefile.am for scuba
#

AM_CPPFLAGS=-I/usr/include/GL $(GLUT_CFLAGS) $(TCL_CFLAGS) $(tcl_includes) -I$(top_srcdir)/include -I$(top_srcdir)/dicom
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
	TclProgressDisplayManager.cpp \
	togl.c

bin_PROGRAMS = scuba
scuba_SOURCES=$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) $(GUI_SRCS) $(LAYER_SRCS) scuba.cpp 
LDADD= $(addprefix $(top_srcdir)/,$(LIBS_MGH)) $(LIBS_TCL_OPENGL)
scuba_LDFLAGS=

foodir=$(prefix)/lib/tcl
foo_SCRIPTS=scuba.tcl tkcon.tcl

check_PROGRAMS=test_ScubaView 
#test_ScubaFrame
# test_ToglManager \
#        test_ScubaLayers \
#        test_DataCollection test_ScubaDataCollections \
#        test_VolumeCollection test_SurfaceCollection \
#        test_TclCommandManager test_PreferencesManager test_IDTracker \
#        test_ScubaTransform test_Volume3 test_ScubaColorLUT \
#        test_ScubaROI test_ScubaROIVolume

TESTS=$(check_PROGRAMS) 

test_ScubaView_SOURCES=test_ScubaView.cpp $(UTILITIES_SRCS) $(BASIC_SRCS) $(GUI_SRCS) $(FSDATA_SRCS)
test_ScubaFrame_SOURCES=$(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) $(GUI_SRCS) $(LAYER_SRCS) test_ScubaFrame.cpp 


#test_ToglManager_SOURCES=test_ToglManager.cpp $(UTILITIES_SRCS) $(BASIC_SRCS) $(GUI_SRCS) $(FSDATA_SRCS)
#test_ScubaLayers_SOURCES=test_ScubaLayers.cpp $(UTILITIES_SRCS) $(BASIC_SRCS) $(GUI_SRCS) $(LAYER_SRCS)
# test_DataCollection_SOURCES=test_DataCollection.cpp $(UTILITIES_SRCS)	\
# $(BASIC_SRCS) $(FSDATA_SRCS)
# test_ScubaDataCollections_SOURCES=test_ScubaDataCollections.cpp $(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) 
# test_VolumeCollection_SOURCES=test_VolumeCollection.cpp $(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDTAT_SRCS) 
# test_SurfaceCollection_SOURCES= test_SurfaceCollection.cpp $(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) 
# test_TclCommandManager_SOURCES=test_TclCommandManager.cpp $(UTILITIES_SRCS) $(BASIC_SRCS) 
# test_PreferencesManager_SOURCES=test_PreferenceManager.cpp $(UTILITIES_SRCS) $(BASIC_SRCS) 
# test_IDTracker_SOURCES=test_IDTracker.cpp
# test_ScubaTransform_SOURCES=test_ScubaTransform.cpp $(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDTAT_SRCS) 
# test_Volume3_SOURCES=test_Volume3.cpp $(UTILITIES_SRCS) $(BASIC_SRCS)
# test_ScubaColorLUT_SOURCES=test_ScubaColorLUT.cpp $(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDTAT_SRCS) 
# test_ScubaROI_SOURCES=test_ScubaROI.cpp $(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDATA_SRCS) 
# test_ScubaROIVolume_SOURCES=test_ScubaROIVolume.cpp $(UTILITIES_SRCS) $(BASIC_SRCS) $(FSDTAT_SRCS) 


