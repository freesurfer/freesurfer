#
# Makefile.am for libsvm.a
#

AM_CPPFLAGS= -Wno-deprecated

AM_LDFLAGS=

noinst_LIBRARIES= libsvm.a

libsvm_a_SOURCES=svm-data-param.cpp svm-kernel.cpp svm-lib-c.cpp \
svm-lib.cpp svm-model.cpp svm-param.cpp svm-similarity.cpp svm-sketch.cpp

am_libqdecprojet_a_OBJECTS=$(libsvm_a_SOURCES:%.cpp=%.o)

EXTRA_DIST=$(wildcard *.h)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
