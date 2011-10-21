#
# Makefile.am
#

noinst_LIBRARIES=libsbl.a

libsbl_a_CXXFLAGS=-I$(top_srcdir)/histo_register/SimpleBaseLib/include/ \
	$(OPENCV_CXXFLAGS) 
libsbl_a_CPPFLAGS=-DUSE_OPENCV -DUSE_ZLIB -DUSE_CDT
libsbl_a_SOURCES=\
	src/system/Timer.cc \
	src/system/FileSystem.cc \
	src/system/SerialPort.cc \
	src/system/TimeUtil.cc \
	src/system/Microcontroller.cc \
	src/system/Signal.cc \
	src/system/Socket.cc \
	src/core/StringUtil.cc \
	src/core/String.cc \
	src/core/Init.cc \
	src/core/Table.cc \
	src/core/Command.cc \
	src/core/Config.cc \
	src/core/UnitTest.cc \
	src/core/File.cc \
	src/core/Display.cc \
	src/core/PathConfig.cc \
	src/other/DrawingLayer.cc \
	src/other/ImageDrawingLayer.cc \
	src/other/SVG.cc \
	src/other/Plot.cc \
	src/other/CodeCheck.cc \
	src/other/Scripting.cc \
	src/other/TaggedFile.cc \
	src/math/ConfigOptimizer.cc \
	src/math/VectorUtil.cc \
	src/math/TensorUtil.cc \
	src/math/Optimizer.cc \
	src/math/MathUtil.cc \
	src/math/Triangulation.cc \
	src/math/KMeans.cc \
	src/math/Geometry.cc \
	src/math/MatrixUtil.cc \
	src/math/TimeSeries.cc \
	src/math/OptimizerUtil.cc \
	src/image/Video.cc \
	src/image/MotionFieldSeq.cc \
	src/image/Track.cc \
	src/image/ImageDraw.cc \
	src/image/ImageRegister.cc \
	src/image/MotionField.cc \
	src/image/Filter.cc \
	src/image/MotionFieldUtil.cc \
	src/image/ImageTransform.cc \
	src/image/Image.cc \
	src/image/ImageSeqUtil.cc \
	src/image/ImageUtil.cc \
  external/CDT/CDT.cc

EXTRA_DIST=\
	external/CDT/CDT.h \
	external/win_dirent.h \
	external/win_timegm.h \
	include/sbl/gui/ConfigEditor.h \
	include/sbl/gui/ImageViewer.h \
	include/sbl/gui/MiscWidgets.h \
	include/sbl/gui/GraphicViewer.h \
	include/sbl/gui/ImageSeqViewer.h \
	include/sbl/gui/MainWindow.h \
	include/sbl/gui/ConfigViewer.h \
	include/sbl/core/Config.h \
	include/sbl/core/StringUtil.h \
	include/sbl/core/Table.h \
	include/sbl/core/Pointer.h \
	include/sbl/core/Init.h \
	include/sbl/core/File.h \
	include/sbl/core/PathConfig.h \
	include/sbl/core/Command.h \
	include/sbl/core/Display.h \
	include/sbl/core/Array.h \
	include/sbl/core/String.h \
	include/sbl/core/Dict.h \
	include/sbl/core/UnitTest.h \
	include/sbl/system/SerialPort.h \
	include/sbl/system/TimeUtil.h \
	include/sbl/system/FileSystem.h \
	include/sbl/system/Signal.h \
	include/sbl/system/Socket.h \
	include/sbl/system/Timer.h \
	include/sbl/system/Microcontroller.h \
	include/sbl/other/CodeCheck.h \
	include/sbl/other/SVG.h \
	include/sbl/other/Scripting.h \
	include/sbl/other/ImageDrawingLayer.h \
	include/sbl/other/TaggedFile.h \
	include/sbl/other/Plot.h \
	include/sbl/other/DrawingLayer.h \
	include/sbl/math/KMeans.h \
	include/sbl/math/Optimizer.h \
	include/sbl/math/Vector.h \
	include/sbl/math/MatrixUtil.h \
	include/sbl/math/Matrix.h \
	include/sbl/math/MathUtil.h \
	include/sbl/math/OptimizerUtil.h \
	include/sbl/math/Tensor.h \
	include/sbl/math/Triangulation.h \
	include/sbl/math/VectorUtil.h \
	include/sbl/math/TimeSeries.h \
	include/sbl/math/TensorUtil.h \
	include/sbl/math/ConfigOptimizer.h \
	include/sbl/math/Geometry.h \
	include/sbl/image/ImageDraw.h \
	include/sbl/image/Video.h \
	include/sbl/image/MotionFieldSeq.h \
	include/sbl/image/Image.h \
	include/sbl/image/Filter.h \
	include/sbl/image/ImageTransform.h \
	include/sbl/image/ImageSeqUtil.h \
	include/sbl/image/MotionField.h \
	include/sbl/image/ImageUtil.h \
	include/sbl/image/MotionFieldUtil.h \
	include/sbl/image/ImageRegister.h \
	include/sbl/image/Track.h

