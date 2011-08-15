#
# Makefile.am
#

noinst_LIBRARIES=libsbl.a

libsbl_a_CXXFLAGS=-I./include/ $(OPENCV_CXXFLAGS) 
libsbl_a_CPPFLAGS=-Wno-conversion
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
	src/image/ImageUtil.cc


