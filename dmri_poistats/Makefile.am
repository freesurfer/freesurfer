#
# Makefile.am for poistats
#

# include this so that the -DHAVE_CONFIG_H flag isn't included when compiling
DEFS=

#itk is required for building poistats
# add ITK
if HAVE_ITK_LIBS

ITK_THINGS=$(ITK_CFLAGS) -DHAVE_ITK_LIBS

UI_DIR = ui
UI_SRCS = \
  $(UI_DIR)/CommandParser.cxx \
  $(UI_DIR)/FreeSurferExecutable.cxx

DATAMODEL_DIR = datamodel
DATAMODEL_SRCS = \
	$(DATAMODEL_DIR)/PoistatsModel.cxx (DATAMODEL_DIR)/PoistatsModel.h

DATA_DIR = $(DATAMODEL_DIR)/data
DATA_SRCS = \
	$(DATA_DIR)/PoistatsReplica.cxx $(DATA_DIR)/PoistatsReplica.h \
	$(DATA_DIR)/PoistatsReplicas.cxx $(DATA_DIR)/PoistatsReplicas.h

EVENTS_DIR = $(DATAMODEL_DIR)/events
EVENTS_SRCS = \
  $(EVENTS_DIR)/PoistatsEvents.h \
  $(EVENTS_DIR)/CommandUpdate.cxx $(EVENTS_DIR)/CommandUpdate.h

UTILS_DIR = $(DATAMODEL_DIR)/utils
UTILS_SRCS = \
	$(UTILS_DIR)/itkPoistatsFilter.txx $(UTILS_DIR)/itkPoistatsFilter.h

BSPLINES_DIR = $(UTILS_DIR)/bspline
BSPLINES_SRCS = \
	$(BSPLINES_DIR)/itkBSplineKernelFunction.txx \
	  $(BSPLINES_DIR)/itkBSplineKernelFunction.h \
	$(BSPLINES_DIR)/itkBSplineScatteredDataPointSetToImageFilter.txx \
    $(BSPLINES_DIR)/itkBSplineScatteredDataPointSetToImageFilter.h \
	$(BSPLINES_DIR)/itkPointSetToImageFilter.txx \
	  $(BSPLINES_DIR)/itkPointSetToImageFilter.h

IO_DIR = io
IO_SRCS = \
	$(IO_DIR)/AsymmetricTensorReaderStrategy.cxx \
	  $(IO_DIR)/AsymmetricTensorReaderStrategy.h \
	$(IO_DIR)/AsymmetricTensorVectorReaderStrategy.cxx \
	  $(IO_DIR)/AsymmetricTensorVectorReaderStrategy.h \
	$(IO_DIR)/SymmetricTensorReaderStrategy.cxx \
	  $(IO_DIR)/SymmetricTensorReaderStrategy.h \
	$(IO_DIR)/TensorReaderStrategy.cxx $(IO_DIR)/TensorReaderStrategy.h \
	$(IO_DIR)/itkMGHImageIO.cxx $(IO_DIR)/itkMGHImageIO.h \
	$(IO_DIR)/itkMGHImageIOFactory.cxx $(IO_DIR)/itkMGHImageIOFactory.h

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
	-lITKNumerics -litkvnl_inst -litkvnl_algo -litkv3p_netlib -litkvnl -litkvcl \
	-lITKMetaIO -litksys -lITKEXPAT -lITKniftiio -lITKznz -litkzlib
dmri_poistats_LDFLAGS = $(OS_LDFLAGS)
dmri_poistats_CXXFLAGS = 

endif
#end itk

EXTRA_DIST=datamodel io

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra

