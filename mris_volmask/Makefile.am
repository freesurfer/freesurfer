##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include $(VTK_CXXFLAGS) -Wno-deprecated

bin_PROGRAMS=mris_volmask mris_volmask_novtk

foodir=$(prefix)/docs/xml
foo_DATA=mris_volmask.help.xml

if HAVE_VTK_LIBS
bin_PROGRAMS += mris_volmask_vtk
mris_volmask_vtk_SOURCES=mris_volmask_old.cpp cmd_line_interface.cpp
mris_volmask_vtk_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH)) \
$(X_LIBS) \
-ldl -lpthread \
$(VTK_LIBS) \
$(VTK_FREETYPE_LIB) $(VTK_PNG_LIB) \
$(VTK_TIFF_LIB) $(VTK_JPEG_LIB) $(VTK_SQLITE_LIB) \
$(VTK_Z_LIB) $(VTK_EXPAT_LIB) $(VTK_METAIO_LIB) \
-lvtkFiltering $(VTK_VERDICT_LIB) -lvtkHybrid $(VTK_HDF5_LIB) \
-lvtkRendering -lvtkGraphics \
-lvtkImaging -lvtkftgl \
-lvtkIO	$(VTK_LSDYNA_LIB) -lvtkFiltering -lvtkCommon -lvtkDICOMParser \
-lvtksys \
$(VTK_MPEG2ENCODE_LIB) -lvtkexoIIc -lvtkNetCDF $(VTK_NETCDFCXX_LIB)
endif

mris_volmask_novtk_SOURCES=mris_volmask_old.cpp cmd_line_interface.cpp
mris_volmask_novtk_CPPFLAGS=-DNO_VTK
mris_volmask_novtk_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))

mris_volmask_SOURCES=mris_volmask.cpp cmd_line_interface.cpp \
	mris_volmask.help.xml.h
mris_volmask_CPPFLAGS=-DNO_VTK
mris_volmask_CXXFLAGS=$(AM_CFLAGS) $(VXL_CFLAGS) -Wno-non-virtual-dtor
mris_volmask_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH)) 

# trick to get test data into the build directory
foo_SOURCES=
foo:
	if ! test -f "$(top_builddir)/mri_aparc2aseg/testdata2.tar.gz"; then \
		cp -v $(top_srcdir)/mri_aparc2aseg/testdata2.tar.gz \
					$(top_builddir)/mri_aparc2aseg/; fi;
#	if ! test -f "$(top_builddir)/scripts/vtk_setup"; then \
#		cp -v $(top_srcdir)/scripts/*_setup \
#					$(top_builddir)/scripts/; fi;
check_PROGRAMS=foo
TESTS=test_mris_volmask $(top_builddir)/scripts/help_xml_validate

# CPPUNIT tests
if HAVE_CPPUNIT

#check_PROGRAMS = TestMRISOBBTree TestMRISDistanceField 

#TESTS = TestMRISOBBTree TestMRISDistanceField

#TestMRISOBBTree_SOURCES = TestMRISOBBTree.cpp  
#TestMRISOBBTree_CXXFLAGS=$(AM_CFLAGS) $(CPPUNIT_CFLAGS) $(VXL_CFLAGS) -Wno-non-virtual-dtor
#TestMRISOBBTree_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(CPPUNIT_LIBS) -lcppunit

#TestMRISDistanceField_SOURCES = TestMRISDistanceField.cpp  
#TestMRISDistanceField_CXXFLAGS=$(AM_CFLAGS) $(CPPUNIT_CFLAGS) $(VXL_CFLAGS) -Wno-non-virtual-dtor
#TestMRISDistanceField_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(CPPUNIT_LIBS) -lcppunit
endif

# put a wrapper around the bin, used to setup VTK enviro vars
if HAVE_VTK_LIBS
install-exec-hook:
	cp $(top_builddir)/mris_volmask/mris_volmask_vtk \
		$(DESTDIR)$(bindir)/mris_volmask_vtk.bin
	echo "#!/bin/tcsh -ef" > $(DESTDIR)$(bindir)/mris_volmask_vtk
	echo "source \$$FREESURFER_HOME/bin/vtk_setup" \
		>> $(DESTDIR)$(bindir)/mris_volmask_vtk
	echo "mris_volmask_vtk.bin \$$argv" >> $(DESTDIR)$(bindir)/mris_volmask_vtk
	chmod a+x $(DESTDIR)$(bindir)/mris_volmask_vtk

uninstall-hook:
	rm -f $(DESTDIR)$(bindir)/mris_volmask_vtk.bin
endif

EXTRA_DIST=cmd_line_interface.h test_mris_volmask $(foo_DATA) \
	TestMRISOBBTree.cpp TestMRISDistanceField.cpp $(BUILT_SOURCES)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra

BUILT_SOURCES=mris_volmask.help.xml.h

