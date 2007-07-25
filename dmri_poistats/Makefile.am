#
# Makefile.am for poistats
#

SUBDIRS=test

#itk is required for building poistats
# add ITK
if HAVE_ITK_LIBS

ITK_THINGS=$(ITK_CFLAGS) -DHAVE_ITK_LIBS

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
	datamodel/utils/itkPoistatsFilter.h \
	datamodel/utils/InitializePath.cxx \
	datamodel/utils/InitializePath.h \
	datamodel/utils/EigenVectorInitPathStrategy.cxx \
	datamodel/utils/EigenVectorInitPathStrategy.h \
	datamodel/utils/FieldLineInitPathStrategy.cxx \
	datamodel/utils/FieldLineInitPathStrategy.h	

IO_SRCS = \
	io/AsymmetricTensorReaderStrategy.cxx \
	io/AsymmetricTensorReaderStrategy.h \
	io/AsymmetricTensorVectorReaderStrategy.cxx \
	io/AsymmetricTensorVectorReaderStrategy.h \
	io/SymmetricTensorReaderStrategy.cxx \
	io/SymmetricTensorReaderStrategy.h \
	io/TensorReaderStrategy.cxx \
	io/TensorReaderStrategy.h

AM_CPPFLAGS =	$(X_CFLAGS) $(VXL_CFLAGS) \
  -I$(top_srcdir)/include \
	-I$(top_srcdir)/itkio \
	-I$(top_srcdir)/itkutils \
  $(ITK_THINGS) -DHAVE_ZLIB

bin_PROGRAMS = \
  dmri_poistats \
  dmri_initialize_path \
  dmri_mask

# poistats
dmri_poistats_SOURCES = \
  $(DATAMODEL_SRCS) \
  $(DATA_SRCS) \
  $(EVENTS_SRCS) \
  $(UTILS_SRCS) \
  $(IO_SRCS) \
  dmri_poistats.cxx 
dmri_poistats_LDADD = $(LIBS_ITK) \
	$(top_srcdir)/itkio/libitkio.a \
	$(top_srcdir)/itkutils/libitkutils.a \
	-lITKBasicFilters -litkgdcm -litkpng -litktiff -lITKSpatialObject \
	-lITKNumerics -litkvnl_inst -litkvnl_algo -litkv3p_netlib \
	-litkvnl -litkvcl \
	-lITKMetaIO -litksys -lITKEXPAT -lITKniftiio -lITKznz -litkzlib \
	$(addprefix $(top_srcdir)/, $(LIBS_MGH) )
dmri_poistats_LDFLAGS = $(OS_LDFLAGS)
dmri_poistats_CXXFLAGS = 

# initialize poistats path stuff
dmri_initialize_path_SOURCES = \
  $(DATAMODEL_SRCS) \
  $(UTILS_SRCS) \
  dmri_initialize_path.cxx 
dmri_initialize_path_LDADD = $(LIBS_ITK) \
	$(top_srcdir)/itkutils/libitkutils.a \
	-lITKBasicFilters -litkgdcm -litkpng -litktiff -lITKSpatialObject \
	-lITKNumerics -litkvnl_inst -litkvnl_algo -litkv3p_netlib \
	-litkvnl -litkvcl \
	-lITKMetaIO -litksys -lITKEXPAT -lITKniftiio -lITKznz -litkzlib \
	$(addprefix $(top_srcdir)/, $(LIBS_MGH) )
dmri_initialize_path_LDFLAGS = $(OS_LDFLAGS)
dmri_initialize_path_CXXFLAGS = 

# mask path
dmri_mask_SOURCES = \
  dmri_mask.cxx 
dmri_mask_LDADD = $(addprefix $(top_srcdir)/, $(LIBS_MGH) ) \
	$(top_srcdir)/itkutils/libitkutils.a
dmri_mask_LDFLAGS = $(OS_LDFLAGS)
dmri_mask_CXXFLAGS = 

endif
#end itk

EXTRA_DIST=datamodel io $(wildcard *.h) $(wildcard *.txx)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra

