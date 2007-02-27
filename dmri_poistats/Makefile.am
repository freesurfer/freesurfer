#
# Makefile.am for poistats
#

#itk is required for building poistats
# add ITK
if HAVE_ITK_LIBS

ITK_THINGS=$(ITK_CFLAGS) -DHAVE_ITK_LIBS

UI_SRCS = \
  ui/CommandParser.cxx \
  ui/FreeSurferExecutable.cxx

DATAMODEL_SRCS = \
	datamodel/PoistatsModel.cxx \
	datamodel/PoistatsModel.h

DATA_SRCS = \
	datamodel/data/PoistatsReplica.cxx \
	datamodel/data/PoistatsReplica.h \
	datamodel/data/PoistatsReplicas.cxx \
	datamodel/data/PoistatsReplicas.h

EVENTS_SRCS = \
  datamodel/events/PoistatsEvents.h \
  datamodel/events/CommandUpdate.cxx \
  datamodel/events/CommandUpdate.h

UTILS_SRCS = \
	datamodel/utils/itkPoistatsFilter.txx \
	datamodel/utils/itkPoistatsFilter.h

BSPLINES_SRCS = \
	datamodel/utils/bspline/itkBSplineKernelFunction.txx \
	datamodel/utils/bspline/itkBSplineKernelFunction.h \
	datamodel/utils/bspline/itkBSplineDataPointSetToImageFilter.txx \
  datamodel/utils/bspline/itkBSplineDataPointSetToImageFilter.h \
	datamodel/utils/bspline/itkPointSetToImageFilter.txx \
	datamodel/utils/bspline/itkPointSetToImageFilter.h

IO_SRCS = \
	io/AsymmetricTensorReaderStrategy.cxx \
	io/AsymmetricTensorReaderStrategy.h \
	io/AsymmetricTensorVectorReaderStrategy.cxx \
	io/AsymmetricTensorVectorReaderStrategy.h \
	io/SymmetricTensorReaderStrategy.cxx \
	io/SymmetricTensorReaderStrategy.h \
	io/TensorReaderStrategy.cxx \
	io/TensorReaderStrategy.h \
	io/itkMGHImageIO.cxx \
	io/itkMGHImageIO.h \
	io/itkMGHImageIOFactory.cxx \
	io/itkMGHImageIOFactory.h

AM_CPPFLAGS =	$(X_CFLAGS) $(VXL_CFLAGS) \
  -I$(top_srcdir)/include \
  $(ITK_THINGS) -DHAVE_ZLIB

bin_PROGRAMS = dmri_poistats
dmri_poistats_SOURCES = \
  $(UI_SRCS) \
  $(BSPLINES_SRCS) \
  $(DATAMODEL_SRCS) \
  $(DATA_SRCS) \
  $(EVENTS_SRCS) \
  $(UTILS_SRCS) \
  $(IO_SRCS) \
  dmri_poistats.cxx 
dmri_poistats_LDADD = $(LIBS_ITK) \
	-lITKBasicFilters -litkgdcm -litkpng -litktiff -lITKSpatialObject \
	-lITKNumerics -litkvnl_inst -litkvnl_algo -litkv3p_netlib \
	-litkvnl -litkvcl \
	-lITKMetaIO -litksys -lITKEXPAT -lITKniftiio -lITKznz -litkzlib
dmri_poistats_LDFLAGS = $(OS_LDFLAGS)
dmri_poistats_CXXFLAGS = 

endif
#end itk

EXTRA_DIST=datamodel io ui $(wildcard *.h) $(wildcard *.txx)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=io/CVS ui/CVS datamodel/CVS datamodel/utils/CVS \
	datamodel/utils/bspline/CVS datamodel/events/CVS datamodel/data/CVS
include $(top_srcdir)/Makefile.extra

