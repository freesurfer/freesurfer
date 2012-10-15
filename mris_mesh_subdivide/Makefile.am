
AM_CFLAGS=-I$(top_srcdir)/include
AM_CPPFLAGS=-I$(top_srcdir)/include \
	-I/usr/include/libxml2
AM_CXXFLAGS=-Wno-deprecated -I$(top_srcdir)/include \
	-I$(top_srcdir)/vtkfsio \
	-I$(top_srcdir)/vtkutils \
	$(VTK_CXXFLAGS)

bin_PROGRAMS = mris_mesh_subdivide
mris_mesh_subdivide_SOURCES=mris_mesh_subdivide.cxx
mris_mesh_subdivide_LDADD= \
	$(top_builddir)/vtkfsio/libvtkfsio.a \
	$(top_builddir)/vtkutils/libvtkutils.a \
	$(addprefix $(top_builddir)/, $(LIBS_MGH))  \
	$(VTK_LIBS) $(VTK_VERDICT_LIB) -lvtkGraphics \
	$(VTK_EXPAT_LIB) $(VTK_FREETYPE_LIB) \
	$(VTK_TIFF_LIB) $(VTK_JPEG_LIB) \
	$(VTK_METAIO_LIB) $(VTK_PNG_LIB) $(VTK_Z_LIB) $(VTK_SQLITE_LIB) \
	-lvtkImaging -lvtkFiltering \
	-lvtkCommon -lvtksys \
	-lvtkexoIIc -lvtkNetCDF $(VTK_NETCDFCXX_LIB) \
	-lvtkRendering -lvtkftgl \
	-lvtkHybrid $(VTK_HDF5_LIB) \
	-lvtkIO $(VTK_LSDYNA_LIB) -lvtkDICOMParser $(VTK_MPEG2ENCODE_LIB)
mris_mesh_subdivide_LDFLAGS=$(OS_LDFLAGS)

# put a script wrapper around the binary, to set up VTK environment (from QDEC example)
install-exec-hook:mris_mesh_subdivide
	cp $(top_builddir)/mris_mesh_subdivide/mris_mesh_subdivide $(DESTDIR)$(bindir)/mris_mesh_subdivide.bin
	echo "#!/bin/tcsh -ef" > $(DESTDIR)$(bindir)/mris_mesh_subdivide
	echo "source \$$FREESURFER_HOME/bin/vtk_setup" \
		>> $(DESTDIR)$(bindir)/mris_mesh_subdivide
	echo "if (-e /bin/pwd) then" >> $(DESTDIR)$(bindir)/mris_mesh_subdivide
	echo "  setenv HOME \`cd \$$HOME; /bin/pwd\`" \
		>> $(DESTDIR)$(bindir)/mris_mesh_subdivide
	echo "  setenv SUBJECTS_DIR \`cd \$$SUBJECTS_DIR; /bin/pwd\`" \
		>> $(DESTDIR)$(bindir)/mris_mesh_subdivide
	echo "endif" >> $(DESTDIR)$(bindir)/mris_mesh_subdivide
	echo "mris_mesh_subdivide.bin \$$argv" >> $(DESTDIR)$(bindir)/mris_mesh_subdivide
	chmod a+x $(DESTDIR)$(bindir)/mris_mesh_subdivide
	if [[ -f /usr/pubsw/bin/chrpath ]]; then \
	/usr/pubsw/bin/chrpath -d $(DESTDIR)$(bindir)/mris_mesh_subdivide.bin; fi

uninstall-hook:
	rm -f $(DESTDIR)$(bindir)/mris_mesh_subdivide.bin



# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
