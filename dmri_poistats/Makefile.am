#
# Makefile.am for poistats
#

#itk is required for building poistats
# add ITK
if ! HAVE_MAC_OSX
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
  -I$(top_srcdir)/dmri_poistats \
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
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKBasicFilters -litkgdcm -litkpng -litktiff -lITKSpatialObject \
	-lITKNumerics -litkvnl_inst -litkvnl_algo -litkv3p_netlib \
	-litkvnl -litkvcl \
	-lITKMetaIO -litksys -lITKEXPAT -lITKniftiio -lITKznz -litkzlib \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
dmri_poistats_LDFLAGS = $(OS_LDFLAGS)
dmri_poistats_CXXFLAGS = 

# initialize poistats path stuff
dmri_initialize_path_SOURCES = \
  $(DATAMODEL_SRCS) \
  $(UTILS_SRCS) \
  dmri_initialize_path.cxx 
dmri_initialize_path_LDADD = $(LIBS_ITK) \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKBasicFilters -litkgdcm -litkpng -litktiff -lITKSpatialObject \
	-lITKNumerics -litkvnl_inst -litkvnl_algo -litkv3p_netlib \
	-litkvnl -litkvcl \
	-lITKMetaIO -litksys -lITKEXPAT -lITKniftiio -lITKznz -litkzlib \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
dmri_initialize_path_LDFLAGS = $(OS_LDFLAGS)
dmri_initialize_path_CXXFLAGS = 

# mask path
dmri_mask_SOURCES = \
  dmri_mask.cxx 
dmri_mask_LDADD = $(addprefix $(top_builddir)/, $(LIBS_MGH) ) \
	$(top_builddir)/itkutils/libitkutils.a
dmri_mask_LDFLAGS = $(OS_LDFLAGS)
dmri_mask_CXXFLAGS = 

# trick to get test data into the build directory
foo_SOURCES=
foo:
	if ! test -f "$(top_builddir)/dmri_poistats/test/data/Path.nii.gz"; then \
		mkdir -p $(top_builddir)/dmri_poistats/test/data/; \
		cp -v $(top_srcdir)/dmri_poistats/test/data/*.nii.gz \
					$(top_builddir)/dmri_poistats/test/data/; fi;

# tests
if HAVE_CPPUNIT

check_PROGRAMS = \
	foo \
  TestPoistatsReplica \
  TestPoistatsReplicas \
  TestPoistats

TESTS = \
	TestPoistatsReplica \
  TestPoistatsReplicas \
  TestPoistats

# TestPoistatsReplica
TestPoistatsReplica_SOURCES= \
  $(DATAMODEL_SRCS) \
  $(DATA_SRCS) \
  $(EVENTS_SRCS) \
  $(UTILS_SRCS) \
  TestPoistatsReplica.cxx

TestPoistatsReplica_CXXFLAGS=$(CPPUNIT_CFLAGS) -Wno-non-virtual-dtor

TestPoistatsReplica_LDADD=$(addprefix $(top_builddir)/, $(LIBS_MGH)) \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	$(CPPUNIT_LIBS) -lcppunit -ldl 	\
	-lITKBasicFilters -litkgdcm -litkpng -litktiff -lITKSpatialObject \
	-lITKNumerics -litkvnl_inst -litkvnl_algo -litkv3p_netlib \
	-litkvnl -litkvcl -lITKMetaIO -litksys -lITKEXPAT -lITKniftiio \
	-lITKznz -litkzlib

# TestPoistatsReplicas
TestPoistatsReplicas_SOURCES= \
  $(DATAMODEL_SRCS) \
  $(DATA_SRCS) \
  $(EVENTS_SRCS) \
  $(UTILS_SRCS) \
  TestPoistatsReplicas.cxx

TestPoistatsReplicas_CXXFLAGS=$(CPPUNIT_CFLAGS) -Wno-non-virtual-dtor

TestPoistatsReplicas_LDADD=$(addprefix $(top_builddir)/, $(LIBS_MGH)) \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	$(CPPUNIT_LIBS) -lcppunit -ldl 	\
	-lITKBasicFilters -litkgdcm -litkpng -litktiff -lITKSpatialObject \
	-lITKNumerics -litkvnl_inst -litkvnl_algo -litkv3p_netlib \
	-litkvnl -litkvcl -lITKMetaIO -litksys -lITKEXPAT -lITKniftiio \
	-lITKznz -litkzlib

# TestPoistats
TestPoistats_SOURCES= \
  $(DATAMODEL_SRCS) \
  $(DATA_SRCS) \
  $(EVENTS_SRCS) \
  $(UTILS_SRCS) \
  $(IO_SRCS) \
  TestPoistats.cxx

TestPoistats_CXXFLAGS=$(CPPUNIT_CFLAGS) -Wno-non-virtual-dtor

TestPoistats_LDADD= \
	$(addprefix $(top_builddir)/, $(LIBS_MGH)) \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	$(LIBS_ITK) \
	$(CPPUNIT_LIBS)	-lcppunit \
	-lITKBasicFilters -litkgdcm -litkpng -litktiff -lITKSpatialObject \
	-lITKNumerics -litkvnl_inst -litkvnl_algo -litkv3p_netlib -litkvnl -litkvcl \
	-lITKIO 	-lITKMetaIO -litksys -lITKEXPAT -lITKniftiio -lITKznz -litkzlib

TestPoistats_LDFLAGS = $(OS_LDFLAGS)
endif
#end cppunit

endif
#end itk

EXTRA_DIST=test/data/FinalSeeds.nii.gz \
	test/data/Path.nii.gz \
	test/data/HelixTensors.nii.gz


# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=io/CVS datamodel/CVS datamodel/utils/CVS datamodel/events/CVS \
datamodel/data/CVS test/CVS test/data/CVS

include $(top_srcdir)/Makefile.extra
endif
