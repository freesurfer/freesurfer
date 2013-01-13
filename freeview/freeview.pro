QT       += core gui

TARGET = freeview
TEMPLATE = app

SOURCES += \
    Annotation2D.cpp \
    CommandEdit.cpp \
    Contour2D.cpp \
    Cursor2D.cpp \
    Cursor3D.cpp \
    CursorFactory.cpp \
    BrushProperty.cpp  \
    DialogAbout.cpp \
    DialogCropVolume.cpp \
    DialogGradientFilter.cpp \
    DialogLoadDTI.cpp \
    DialogLoadPointSet.cpp \
    DialogLoadVolume.cpp \
    DialogNewROI.cpp \
    DialogNewPointSet.cpp \
    DialogNewVolume.cpp \
    DialogPreferences.cpp \
    DialogSavePointSet.cpp \
    DialogSaveScreenshot.cpp \
    DialogTransformVolume.cpp \
    DialogVolumeFilter.cpp \
    DialogWriteMovieFrames.cpp \
    FloatingStatusBar.cpp \
    FSLabel.cpp\
    FSPointSet.cpp \
    FSSurface.cpp \
    FSVolume.cpp \
    GenericRenderView.cpp \
    InfoTreeWidget.cpp \
    Interactor.cpp \
    Interactor2D.cpp \
    Interactor2DMeasure.cpp \
    Interactor2DNavigate.cpp \
    Interactor2DPointSetEdit.cpp \
    Interactor2DROIEdit.cpp \
    Interactor2DVolumeCrop.cpp \
    Interactor2DVolumeEdit.cpp \
    Interactor2DVoxelEdit.cpp \
    Interactor3D.cpp \
    Interactor3DMeasure.cpp \
    Interactor3DNavigate.cpp \
    Interactor3DVolumeCrop.cpp \
    Layer.cpp \
    LayerCollection.cpp \
    LayerDTI.cpp \
    LayerEditable.cpp \
    LayerMRI.cpp \
    LayerPLabel.cpp \
    LayerPointSet.cpp \
    LayerProperty.cpp \
    LayerPropertyDTI.cpp \
    LayerPropertyMRI.cpp \
    LayerPropertyPointSet.cpp \
    LayerPropertyROI.cpp \
    LayerPropertySurface.cpp \
    LayerROI.cpp \
    LayerSurface.cpp \
    LayerTreeWidget.cpp \
    LayerVolumeBase.cpp \
    LivewireTool.cpp \
    LUTDataHolder.cpp \
    main.cpp \
    MainWindow.cpp \
    MyCmdLineParser.cpp \
    MyUtils.cpp \
    MyVTKUtils.cpp \
    PanelLayer.cpp \
    PanelPointSet.cpp \
    PanelROI.cpp \
    PanelSurface.cpp \
    PanelVolume.cpp \
    qtcolorpicker.cpp \
    QVTKWidget.cxx \
    QVTKPaintEngine.cxx \
    Region2D.cpp \
    Region2DLine.cpp \
    Region2DPolyline.cpp \
    Region2DRectangle.cpp \
    RenderView.cpp \
    RenderView2D.cpp \
    RenderView3D.cpp \
    SurfaceAnnotation.cpp \
    SurfaceLabel.cpp \
    SurfaceOverlay.cpp \
    SurfaceOverlayProperty.cpp \
    SurfaceRegion.cpp \
    SurfaceRegionGroups.cpp \
    TermWidget.cpp \
    ThreadBuildContour.cpp \
    ThreadIOWorker.cpp \
    ToolWindowEdit.cpp \
    ToolWindowMeasure.cpp \
    ToolWindowROIEdit.cpp \
    UIUpdateHelper.cpp \
    VolumeCropper.cpp \
    VolumeFilter.cpp \
    VolumeFilterConvolve.cpp \
    VolumeFilterGradient.cpp \
    VolumeFilterMean.cpp \
    VolumeFilterMedian.cpp \
    VolumeFilterSobel.cpp \
    vtkSimpleLabelEdgeFilter.cpp \
    WidgetHistogram.cpp \
    WindowConfigureOverlay.cpp \
    WindowQuickReference.cpp \
    FSTrack.cpp \
    track_io/TrackIO.cpp \
    TrackData.cpp \
    Track.cpp \
    LayerTrack.cpp \
    TrackGroup.cpp \
    PanelTrack.cpp \
    LayerPropertyTrack.cpp \
    DialogSaveVolume.cpp \
    DialogReplaceLabel.cpp \
    LayerVolumeTrack.cpp \
    LayerLandmarks.cpp \
    SurfaceROI.cpp \
    ProgressCallback.cpp \
    MainApplication.cpp \
    DialogRepositionSurface.cpp \
    WindowTimeCourse.cpp \
    WidgetTimeCoursePlot.cpp \
    LayerMRIWorkerThread.cpp \
    DialogLabelStats.cpp \
    VolumeFilterWorkerThread.cpp \
    FSGroupDescriptor.cpp \
    WindowGroupPlot.cpp \
    WidgetGroupPlot.cpp \
    SurfaceSpline.cpp \
    DialogLoadSurfaceOverlay.cpp \
    DialogReloadLayer.cpp \
    DialogSmoothSurface.cpp \
    DialogLineProfile.cpp \
    LayerLineProfile.cpp \
    LayerPropertyLineProfile.cpp

HEADERS  += \
    Annotation2D.h \
    BrushProperty.h \
    CommandEdit.h \
    Contour2D.h \
    Cursor2D.h \
    Cursor3D.h \
    CursorFactory.h \
    DialogAbout.h \
    DialogCropVolume.h \
    DialogGradientFilter.h \
    DialogLoadDTI.h \
    DialogLoadPointSet.h \
    DialogLoadVolume.h \
    DialogPreferences.h \
    DialogNewPointSet.h \
    DialogNewROI.h \
    DialogNewVolume.h \
    DialogSavePointSet.h \
    DialogSaveScreenshot.h \
    DialogTransformVolume.h \
    DialogVolumeFilter.h \
    DialogWriteMovieFrames.h \
    FloatingStatusBar.h \
    FSLabel.h\
    FSPointSet.h \
    FSSurface.h \
    FSVolume.h \
    GenericRenderView.h \
    InfoTreeWidget.h \
    Interactor.h \
    Interactor2D.h \
    Interactor2DMeasure.h \
    Interactor2DNavigate.h \
    Interactor2DPointSetEdit.h \
    Interactor2DROIEdit.h \
    Interactor2DVolumeCrop.h \
    Interactor2DVolumeEdit.h \
    Interactor2DVoxelEdit.h \
    Interactor3D.h \
    Interactor3DMeasure.h \
    Interactor3DNavigate.h \
    Interactor3DVolumeCrop.h \
    Layer.h \
    LayerCollection.h \
    LayerDTI.h \
    LayerEditable.h \
    LayerPLabel.h \
    LayerPointSet.h \
    LayerProperty.h \
    LayerPropertyDTI.h \
    LayerPropertyMRI.h \
    LayerPropertyPointSet.h \
    LayerPropertyROI.h \
    LayerPropertySurface.h \
    LayerROI.h \
    LayerSurface.h \
    LayerTreeWidget.h \
    LayerVolumeBase.h \
    LayerMRI.h \
    LivewireTool.h \
    LUTDataHolder.h \
    MainWindow.h \
    MyCmdLineParser.h \
    MyUtils.h \
    PanelLayer.h \
    PanelPointSet.h \
    PanelROI.h \
    PanelSurface.h \
    PanelVolume.h \
    qtcolorpicker.h \
    QVTKWidget.h \
    Region2D.h \
    Region2DLine.h \
    Region2DPolyline.h \
    Region2DRectangle.h \
    RenderView.h \
    RenderView2D.h \
    RenderView3D.h \
    SurfaceLabel.h \
    SurfaceAnnotation.h \
    SurfaceOverlay.h \
    SurfaceOverlayProperty.h \
    SurfaceRegion.h \
    SurfaceRegionGroups.h \
    TermWidget.h \
    ThreadBuildContour.h \
    ThreadIOWorker.h \
    ToolWindowEdit.h \
    ToolWindowMeasure.h \
    ToolWindowROIEdit.h \
    UIUpdateHelper.h \
    VolumeCropper.h \
    VolumeFilter.h \
    VolumeFilterSobel.h \
    vtkSimpleLabelEdgeFilter.h \
    WidgetHistogram.h \
    WindowConfigureOverlay.h \
    WindowQuickReference.h \
    FSTrack.h \
    TrackData.h \
    LayerTrack.h \
    TrackGroup.h \
    PanelTrack.h \
    LayerPropertyTrack.h \
    DialogSaveVolume.h \
    DialogReplaceLabel.h \
    LayerVolumeTrack.h \
    LayerLandmarks.h \
    SurfaceROI.h \
    ProgressCallback.h \
    MainApplication.h \
    DialogRepositionSurface.h \
    WindowTimeCourse.h \
    WidgetTimeCoursePlot.h \
    LayerMRIWorkerThread.h \
    DialogLabelStats.h \
    VolumeFilterWorkerThread.h \
    FSGroupDescriptor.h \
    WindowGroupPlot.h \
    WidgetGroupPlot.h \
    SurfaceSpline.h \
    DialogLoadSurfaceOverlay.h \
    DialogReloadLayer.h \
    DialogSmoothSurface.h \
    DialogLineProfile.h \
    LayerLineProfile.h \
    LayerPropertyLineProfile.h

FORMS    += MainWindow.ui \
    PanelVolume.ui \
    PanelSurface.ui \
    PanelROI.ui \
    PanelPointSet.ui \
    DialogLoadVolume.ui \
    DialogPreferences.ui \
    ToolWindowMeasure.ui \
    ToolWindowEdit.ui \
    DialogNewVolume.ui \
    DialogNewROI.ui \
    DialogNewPointSet.ui \
    ToolWindowROIEdit.ui \
    DialogLoadPointSet.ui \
    DialogLoadDTI.ui \
    WindowConfigureOverlay.ui \
    DialogTransformVolume.ui \
    DialogCropVolume.ui \
    DialogSaveScreenshot.ui \
    DialogVolumeFilter.ui \
    DialogAbout.ui \
    DialogWriteMovieFrames.ui \
    DialogGradientFilter.ui \
    WindowQuickReference.ui \
    FloatingStatusBar.ui \
    TermWidget.ui \
    PanelTrack.ui \
    DialogSavePointSet.ui \
    DialogSaveVolume.ui \
    DialogReplaceLabel.ui \
    DialogRepositionSurface.ui \
    WindowTimeCourse.ui \
    DialogLabelStats.ui \
    WindowGroupPlot.ui \
    DialogLoadSurfaceOverlay.ui \
    DialogReloadLayer.ui \
    DialogSmoothSurface.ui \
    DialogLineProfile.ui

RESOURCES += \
    freeview.qrc

LIBS += \
    -lvtkverdict -lvtkGraphics -lvtkmetaio -lvtkpng -lvtkzlib \
    -lvtksqlite -lvtkImaging -lvtkFiltering -lvtkCommon -lvtksys \
    -lvtkGenericFiltering -lvtkexoIIc -lvtkNetCDF -lvtkVolumeRendering \
    -lvtkRendering -lvtkftgl -lvtkWidgets -lvtkHybrid -lvtkIO -lvtkDICOMParser

#LIBS += \
#    -lvtkhdf5_hl -lvtkhdf5 -lLSDyna  -lvtkNetCDF_cxx

QMAKE_CXXFLAGS += -Wno-deprecated -DUNICODE -D_FILE_OFFSET_BITS=64 -D_LARGE_FILES \
                                  -Wno-write-strings  -DDEVELOPMENT

# set this to your local dev directory
FREESURFER_DEV_DIR = /homes/5/rpwang/freesurfer/dev

# set this to your local install bin directory
# freeview.bin will be copied to that directory
FREESURFER_BIN = /homes/5/rpwang/freesurfer/bin

# for linux
unix {
INCLUDEPATH += /usr/pubsw/packages/vtk/current/include/vtk-5.6 \
               $$FREESURFER_DEV_DIR/include $$FREESURFER_DEV_DIR/vtkutils \
               /usr/pubsw/packages/mni/current/include \
               $$FREESURFER_DEV_DIR/lineprof

QMAKE_CXXFLAGS += -I/usr/pubsw/packages/itk/current/include/InsightToolkit \
    -I/usr/pubsw/packages/itk/current/include/InsightToolkit/Algorithms \
    -I/usr/pubsw/packages/itk/current/include/InsightToolkit/BasicFilters \
    -I/usr/pubsw/packages/itk/current/include/InsightToolkit/Common \
    -I/usr/pubsw/packages/itk/current/include/InsightToolkit/IO \
    -I/usr/pubsw/packages/itk/current/include/InsightToolkit/Numerics \
    -I/usr/pubsw/packages/itk/current/include/InsightToolkit/Numerics/Statistics \
    -I/usr/pubsw/packages/itk/current/include/InsightToolkit/Review \
    -I/usr/pubsw/packages/itk/current/include/InsightToolkit/Review/Statistics \
    -I/usr/pubsw/packages/itk/current/include/InsightToolkit/SpatialObject \
    -I/usr/pubsw/packages/itk/current/include/InsightToolkit/Utilities \
    -I/usr/pubsw/packages/vxl/current/include/vxl/core \
    -I/usr/pubsw/packages/vxl/current/include/vxl/vcl \
    -I/usr/pubsw/packages/vxl/current/include/vxl/v3p/netlib \
    -I/usr/pubsw/packages/vxl/current/include/vxl/v3p/netlib/opt \
    -I/usr/pubsw/packages/petsc/current/include


LIBS += -L/usr/pubsw/packages/vtk/current/lib/vtk-5.6 -L/usr/X11R6/lib \
    -lX11 -lXext -lXt -lSM -lICE -lGLU -lm -ldl \
    -L/usr/pubsw/packages/vxl/current/lib -L/usr/pubsw/packages/itk/current/lib/InsightToolkit \
    $$FREESURFER_DEV_DIR/utils/libutils.a $$FREESURFER_DEV_DIR/fsgdf/libfsgdf.a \
    $$FREESURFER_DEV_DIR/vtkutils/libvtkutils.a \
    $$FREESURFER_DEV_DIR/lineprof/liblineprof.a \
    $$FREESURFER_DEV_DIR/hipsstubs/libhipsstubs.a $$FREESURFER_DEV_DIR/vtkutils/libvtkutils.a \
    $$FREESURFER_DEV_DIR/rgb/librgb.a $$FREESURFER_DEV_DIR/unix/libunix.a $$FREESURFER_DEV_DIR/dicom/libdicom.a \
    $$FREESURFER_DEV_DIR/jpeg/libjpeg.a $$FREESURFER_DEV_DIR/tiff/libtiff.a $$FREESURFER_DEV_DIR/expat/libexpat.a \
    /usr/pubsw/packages/itk/current/lib/InsightToolkit/libITKIO.a /usr/pubsw/packages/itk/current/lib/InsightToolkit/libITKAlgorithms.a \
    /usr/pubsw/packages/itk/current/lib/InsightToolkit/libITKCommon.a /usr/pubsw/packages/itk/current/lib/InsightToolkit/libITKMetaIO.a \
    /usr/pubsw/packages/itk/current/lib/InsightToolkit/libITKniftiio.a /usr/pubsw/packages/itk/current/lib/InsightToolkit/libITKNrrdIO.a \
    /usr/pubsw/packages/itk/current/lib/InsightToolkit/libitkpng.a /usr/pubsw/packages/itk/current/lib/InsightToolkit/libitksys.a \
    /usr/pubsw/packages/itk/current/lib/InsightToolkit/libitktiff.a /usr/pubsw/packages/itk/current/lib/InsightToolkit/libitkv3p_netlib.a \
    /usr/pubsw/packages/itk/current/lib/InsightToolkit/libitkzlib.a /usr/pubsw/packages/itk/current/lib/InsightToolkit/libitkgdcm.a \
    /usr/pubsw/packages/itk/current/lib/InsightToolkit/libitkopenjpeg.a /usr/pubsw/packages/itk/current/lib/InsightToolkit/libitkjpeg8.a \
    /usr/pubsw/packages/itk/current/lib/InsightToolkit/libitkjpeg12.a /usr/pubsw/packages/itk/current/lib/InsightToolkit/libitkjpeg16.a \
    /usr/pubsw/packages/itk/current/lib/InsightToolkit/libITKDICOMParser.a \
    /usr/lib64/libuuid.a -lz -lcrypt -ldl -lpthread \
    /usr/pubsw/packages/mni/1.4/lib/libvolume_io.a -L/usr/pubsw/packages/mni/1.4/lib /usr/pubsw/packages/mni/1.4/lib/libminc.a /usr/pubsw/packages/mni/1.4/lib/libnetcdf.a \
    -lvnl_algo -lvnl -lvcl -lnetlib -lv3p_netlib \
    -L/usr/pubsw/packages/petsc/current/lib -lpetscts -lpetscsnes -lpetscksp \
    -lpetscdm -lpetscmat -lpetscvec -lpetsc -lmpich -lfmpich \
    /usr/lib64/liblapack.a /usr/lib64/libblas.a -lgfortran

TARGET = freeview.bin
DESTDIR = $$FREESURFER_BIN
}

# for mac
macx {

TARGET = FreeView
RC_FILE = resource/icons/freeview.icns

# uncomment following lines to build for 10.5 compatible binaries
CONFIG += i386
QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.5
QMAKE_MAC_SDK=/Developer/SDKs/MacOSX10.5.sdk
QMAKE_CXXFLAGS += -mmacosx-version-min=10.5 -arch i386
QMAKE_CFLAGS += -mmacosx-version-min=10.5 -arch i386
QMAKE_LFLAGS += -mmacosx-version-min=10.5 -arch i386

LIBS -= -L/usr/pubsw/packages/vtk/current/lib/vtk-5.6 -L/usr/X11R6/lib \
    -L/usr/pubsw/packages/vxl/current/lib -L/usr/pubsw/packages/itk/current/lib/InsightToolkit \
    $$FREESURFER_DEV_DIR/utils/libutils.a $$FREESURFER_DEV_DIR/fsgdf/libfsgdf.a \
    $$FREESURFER_DEV_DIR/hipsstubs/libhipsstubs.a $$FREESURFER_DEV_DIR/vtkutils/libvtkutils.a \
    $$FREESURFER_DEV_DIR/rgb/librgb.a $$FREESURFER_DEV_DIR/unix/libunix.a $$FREESURFER_DEV_DIR/dicom/libdicom.a \
    $$FREESURFER_DEV_DIR/jpeg/libjpeg.a $$FREESURFER_DEV_DIR/tiff/libtiff.a $$FREESURFER_DEV_DIR/expat/libexpat.a \
    /usr/pubsw/packages/itk/current/lib/InsightToolkit/libITKIO.a /usr/pubsw/packages/itk/current/lib/InsightToolkit/libITKAlgorithms.a \
    /usr/pubsw/packages/itk/current/lib/InsightToolkit/libITKCommon.a /usr/pubsw/packages/itk/current/lib/InsightToolkit/libITKMetaIO.a \
    /usr/pubsw/packages/itk/current/lib/InsightToolkit/libITKniftiio.a /usr/pubsw/packages/itk/current/lib/InsightToolkit/libITKNrrdIO.a \
    /usr/pubsw/packages/itk/current/lib/InsightToolkit/libitkpng.a /usr/pubsw/packages/itk/current/lib/InsightToolkit/libitksys.a \
    /usr/pubsw/packages/itk/current/lib/InsightToolkit/libitktiff.a /usr/pubsw/packages/itk/current/lib/InsightToolkit/libitkv3p_netlib.a \
    /usr/pubsw/packages/itk/current/lib/InsightToolkit/libitkzlib.a /usr/pubsw/packages/itk/current/lib/InsightToolkit/libitkgdcm.a \
    /usr/pubsw/packages/itk/current/lib/InsightToolkit/libitkopenjpeg.a /usr/pubsw/packages/itk/current/lib/InsightToolkit/libitkjpeg8.a \
    /usr/pubsw/packages/itk/current/lib/InsightToolkit/libitkjpeg12.a /usr/pubsw/packages/itk/current/lib/InsightToolkit/libitkjpeg16.a \
    /usr/pubsw/packages/itk/current/lib/InsightToolkit/libITKDICOMParser.a /usr/lib64/libuuid.a -lz -lcrypt -ldl -lpthread \
    /usr/pubsw/packages/mni/1.4/lib/libvolume_io.a -L/usr/pubsw/packages/mni/1.4/lib /usr/pubsw/packages/mni/1.4/lib/libminc.a /usr/pubsw/packages/mni/1.4/lib/libnetcdf.a \
    -lvnl_algo -lvnl -lvcl -lnetlib -lv3p_netlib

LIBS -= \
    -lvtkverdict -lvtkGraphics -lvtkmetaio -lvtkpng \ #-lvtkzlib \
    -lvtksqlite -lvtkImaging -lvtkFiltering -lvtkCommon -lvtksys \
    -lvtkGenericFiltering -lvtkexoIIc -lvtkNetCDF -lvtkVolumeRendering \
    -lvtkRendering -lvtkftgl -lvtkWidgets -lvtkHybrid -lvtkIO -lvtkDICOMParser

LIBS += -lvtkHybrid -lvtkVolumeRendering -lvtkRendering -lvtkIO \
          -lvtkGraphics -lvtkGenericFiltering  \
          -lvtkImaging -lvtkDICOMParser -lvtkFiltering -lvtktiff \
          -lvtkCommon -lvtkftgl -lvtkfreetype -lvtkexpat -lvtkjpeg -lvtkpng -lvtksys \
          -lvtkzlib

FREESURFER_DEV_DIR = /Users/rpwang/freesurfer/dev

INCLUDEPATH += /usr/pubsw/packages/vtk/current/include/vtk-5.6 $$FREESURFER_DEV_DIR/include $$FREESURFER_DEV_DIR/vtkutils \
               "/usr/pubsw/packages/mni/current/include"

LIBS += -L/usr/pubsw/packages/vtk/current/lib/vtk-5.6 -framework OpenGL -lm -ldl -lz -framework ApplicationServices \
    -framework CoreServices -framework cocoa -framework IOKit \
    -L/usr/pubsw/packages/vxl/current/lib -L/usr/pubsw/packages/itk/current/lib/InsightToolkit \
    $$FREESURFER_DEV_DIR/lib/libutils.a $$FREESURFER_DEV_DIR/lib/libfsgdf.a \
    $$FREESURFER_DEV_DIR/lib/libhipsstubs.a $$FREESURFER_DEV_DIR/lib/libvtkutils.a \
    $$FREESURFER_DEV_DIR/lib/librgb.a $$FREESURFER_DEV_DIR/lib/libunix.a $$FREESURFER_DEV_DIR/lib/libdicom.a \
    $$FREESURFER_DEV_DIR/lib/libjpeg.a $$FREESURFER_DEV_DIR/lib/libtiff.a $$FREESURFER_DEV_DIR/lib/libexpat.a \
    /usr/pubsw/packages/itk/current/lib/InsightToolkit/libITKIO.a /usr/pubsw/packages/itk/current/lib/InsightToolkit/libITKAlgorithms.a \
    /usr/pubsw/packages/itk/current/lib/InsightToolkit/libITKCommon.a /usr/pubsw/packages/itk/current/lib/InsightToolkit/libITKMetaIO.a \
    /usr/pubsw/packages/itk/current/lib/InsightToolkit/libITKniftiio.a /usr/pubsw/packages/itk/current/lib/InsightToolkit/libITKNrrdIO.a \
    /usr/pubsw/packages/itk/current/lib/InsightToolkit/libitkpng.a /usr/pubsw/packages/itk/current/lib/InsightToolkit/libitksys.a \
    /usr/pubsw/packages/itk/current/lib/InsightToolkit/libitktiff.a /usr/pubsw/packages/itk/current/lib/InsightToolkit/libitkv3p_netlib.a \
    /usr/pubsw/packages/itk/current/lib/InsightToolkit/libitkzlib.a /usr/pubsw/packages/itk/current/lib/InsightToolkit/libitkgdcm.a \
    /usr/pubsw/packages/itk/current/lib/InsightToolkit/libitkopenjpeg.a /usr/pubsw/packages/itk/current/lib/InsightToolkit/libitkjpeg8.a \
    /usr/pubsw/packages/itk/current/lib/InsightToolkit/libitkjpeg12.a /usr/pubsw/packages/itk/current/lib/InsightToolkit/libitkjpeg16.a \
    /usr/pubsw/packages/itk/current/lib/InsightToolkit/libITKDICOMParser.a -lz -ldl -lpthread \
    /usr/pubsw/packages/mni/current/lib/libvolume_io.a -L/usr/pubsw/packages/mni/current/lib \
    /usr/pubsw/packages/mni/current/lib/libminc.a /usr/pubsw/packages/mni/current/lib/libnetcdf.a \
    -lvnl_algo -lvnl -lvcl -lnetlib -lv3p_netlib

LIBS -= -L/usr/X11R6/lib -lX11 -lXext -lXt -lSM -lICE -lGLU -lGL

DESTDIR = ./
}

OTHER_FILES += \
    resource/QuickRef.html \
    Makefile.am
