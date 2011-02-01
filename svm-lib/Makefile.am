#
# Makefile.am for libsvm.a
#

AM_CXXFLAGS= -Wno-deprecated -Wno-unused

AM_LDFLAGS=

noinst_LIBRARIES= libsvm.a

libsvm_a_SOURCES=svm-data-param.cpp svm-kernel.cpp svm-lib-c.cpp \
	svm-lib.cpp svm-model.cpp svm-param.cpp svm-similarity.cpp \
	svm-sketch.cpp

am_libqdecprojet_a_OBJECTS=$(libsvm_a_SOURCES:%.cpp=%.o)

EXTRA_DIST=svm-data-param.h svm-element-types.h svm-kernel-param.h \
	svm-lib-c.h  svm-matrix.h  svm-param.h \
	svm-similarity.h svm-vector-types.h \
	svm-element-type.h svm-io-format.h svm-kernel.h \
	svm-lib.h svm-model.h svm-pr-loqo.h svm-sketch.h \
	svm-vector.h svm-pr-loqo.c \
	classification_setup.example INSTALL_README.txt README.txt

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
