#-------------------------------------------------
#
# Sample qt project file.
# This project is not needed to build within freesurfer
# environment. But it is *strongly* recommended for local
# developing using QtCreator.
#
#-------------------------------------------------

# set this to your local dev directory
FREESURFER_DEV_DIR = /homes/5/rpwang/freesurfer/dev

# set this to your local install bin directory
# freeview.bin will be copied to that directory
FREESURFER_BIN = /homes/5/rpwang/freesurfer/bin

QT       += core gui

TARGET = dummy_qt
TEMPLATE = app


SOURCES += main.cpp\
        MainWindow.cpp

HEADERS  += MainWindow.h

FORMS    += MainWindow.ui

# for linux
unix {

INCLUDEPATH += /usr/pubsw/packages/mni/current/include

LIBS += -L/usr/X11R6/lib \
    -lX11 -lXext -lXt -lSM -lICE -lGLU -ldl \
    -L/usr/pubsw/packages/vxl/current/lib -L/usr/pubsw/packages/itk/current/lib/InsightToolkit \
    $$FREESURFER_DEV_DIR/utils/libutils.a $$FREESURFER_DEV_DIR/fsgdf/libfsgdf.a \
    $$FREESURFER_DEV_DIR/hipsstubs/libhipsstubs.a \
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
    -lvnl_algo -lvnl -lvcl -lnetlib -lv3p_netlib

TARGET = dummy_qt.bin
DESTDIR = $$FREESURFER_BIN
}

# for mac
macx {

#RC_FILE = resource/icons/freeview.icns

# uncomment following lines to build for 10.5 compatible binaries
CONFIG += i386
QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.5
QMAKE_MAC_SDK=/Developer/SDKs/MacOSX10.5.sdk
QMAKE_CXXFLAGS += -mmacosx-version-min=10.5 -arch i386
QMAKE_CFLAGS += -mmacosx-version-min=10.5 -arch i386
QMAKE_LFLAGS += -mmacosx-version-min=10.5 -arch i386

LIBS -= -L/usr/X11R6/lib \
    -L/usr/pubsw/packages/vxl/current/lib -L/usr/pubsw/packages/itk/current/lib/InsightToolkit \
    $$FREESURFER_DEV_DIR/utils/libutils.a $$FREESURFER_DEV_DIR/fsgdf/libfsgdf.a \
    $$FREESURFER_DEV_DIR/hipsstubs/libhipsstubs.a \
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


FREESURFER_DEV_DIR = /Users/rpwang/freesurfer/dev

INCLUDEPATH += "/usr/pubsw/packages/mni/current/include"

LIBS += -framework OpenGL -ldl -lz -framework ApplicationServices \
    -framework CoreServices -framework cocoa -framework IOKit \
    -L/usr/pubsw/packages/vxl/current/lib -L/usr/pubsw/packages/itk/current/lib/InsightToolkit \
    $$FREESURFER_DEV_DIR/lib/libutils.a $$FREESURFER_DEV_DIR/lib/libfsgdf.a \
    $$FREESURFER_DEV_DIR/lib/libhipsstubs.a \
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
    Makefile.am

RESOURCES += \
    dummy_qt.qrc
