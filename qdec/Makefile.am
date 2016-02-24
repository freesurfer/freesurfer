foodir=$(prefix)/lib/images

foo_datafiles=qdec_splash.png

AM_CPPFLAGS=-I$(top_srcdir)/include $(TCL_CFLAGS) \
	-I/usr/include/libxml2

AM_CXXFLAGS=-Wno-deprecated -I$(top_srcdir)/include \
	-I$(top_srcdir)/vtkfsio \
	-I$(top_srcdir)/vtkutils \
	-I$(top_srcdir)/qdecproject \
	$(VTK_CXXFLAGS) $(KWWIDGETS_CXXFLAGS) $(TCL_CFLAGS)

dist_bin_SCRIPTS=run_mris_preproc

if ENABLE_TCLTK_APPS
bin_PROGRAMS = qdec
qdec_SOURCES= \
	QdecMain.cxx \
	vtkKWQdecApp.cxx vtkKWQdecApp.h vtkKWQdecAppTcl.cxx \
	vtkKWQdecWindow.cxx vtkKWQdecWindow.h vtkKWQdecWindowTcl.cxx \
	vtkKWQdecView.cxx vtkKWQdecView.h vtkKWQdecViewTcl.cxx \
	QdecUtilities.cxx QdecLibInit.cxx \
	QdecEvents.h QdecVertexAnnotationLookup.h
qdec_LDADD=  \
	$(top_builddir)/fsgdf/libfsgdf.a \
	$(top_builddir)/vtkfsio/libvtkfsio.a \
	$(top_builddir)/vtkutils/libvtkutils.a \
	$(top_builddir)/qdecproject/libqdecproject.a \
	$(addprefix $(top_builddir)/, $(LIBS_MGH)) $(LIBS_TCL_OPENGL) \
	$(VTK_LIBS) $(VTK_VERDICT_LIB) -lvtkGraphics \
	$(VTK_EXPAT_LIB) $(VTK_FREETYPE_LIB) \
	$(VTK_TIFF_LIB) $(VTK_JPEG_LIB) \
	$(VTK_METAIO_LIB) $(VTK_PNG_LIB) $(VTK_Z_LIB) \
	$(VTK_SQLITE_LIB) \
	-lvtkImaging -lvtkFiltering \
	-lvtkCommon -lvtksys \
	-lvtkexoIIc -lvtkNetCDF $(VTK_NETCDFCXX_LIB) \
	-lvtkVolumeRendering -lvtkRendering -lvtkftgl \
	-lvtkWidgets -lvtkHybrid $(VTK_HDF5_LIB) \
	-lvtkIO	$(VTK_LSDYNA_LIB) -lvtkDICOMParser $(VTK_MPEG2ENCODE_LIB) \
	-lvtkImagingTCL -lvtkVolumeRenderingTCL -lvtkRenderingTCL \
	-lvtkFilteringTCL -lvtkWidgetsTCL \
	-lvtkHybridTCL -lvtkGraphicsTCL -lvtkImagingTCL -lvtkIOTCL \
	-lvtkCommonTCL \
	$(KWWIDGETS_LIBS) -lKWWidgets
qdec_LDFLAGS=$(OS_LDFLAGS) $(LDFLAGS_TCL)

qdecdir=$(prefix)/lib/resource
qdec_DATA=QdecIcons.txt

vtkKWQdecAppTcl.cxx: vtkKWQdecApp.h
	${VTK_DIR}/bin/vtkWrapTcl vtkKWQdecApp.h \
	${VTK_LIB_VTK}/hints 1 vtkKWQdecAppTcl.cxx

vtkKWQdecWindowTcl.cxx: vtkKWQdecWindow.h
	${VTK_DIR}/bin/vtkWrapTcl vtkKWQdecWindow.h \
	${VTK_LIB_VTK}/hints 1 vtkKWQdecWindowTcl.cxx

vtkKWQdecViewTcl.cxx: vtkKWQdecView.h
	${VTK_DIR}/bin/vtkWrapTcl vtkKWQdecView.h \
	${VTK_LIB_VTK}/hints 1 vtkKWQdecViewTcl.cxx

QdecLibInit.cxx: QdecLibInit.data
	${VTK_DIR}/bin/vtkWrapTclInit QdecLibInit.data QdecLibInit.cxx

# put a wrapper around the bin, used to setup tcltktixblt,vtk,kww enviro vars
install-exec-hook:qdec
	cp $(top_builddir)/qdec/qdec $(DESTDIR)$(bindir)/qdec.bin
	echo "#!/bin/tcsh -ef" > $(DESTDIR)$(bindir)/qdec
	echo "source \$$FREESURFER_HOME/bin/tcl_setup" \
		>> $(DESTDIR)$(bindir)/qdec
	echo "source \$$FREESURFER_HOME/bin/vtk_setup" \
		>> $(DESTDIR)$(bindir)/qdec
	echo "source \$$FREESURFER_HOME/bin/kww_setup" \
		>> $(DESTDIR)$(bindir)/qdec
	echo "if (-e /bin/pwd) then" >> $(DESTDIR)$(bindir)/qdec
	echo "  setenv HOME \`cd \$$HOME; /bin/pwd\`" \
		>> $(DESTDIR)$(bindir)/qdec
	echo "  setenv SUBJECTS_DIR \`cd \$$SUBJECTS_DIR; /bin/pwd\`" \
		>> $(DESTDIR)$(bindir)/qdec
	echo "endif" >> $(DESTDIR)$(bindir)/qdec
	echo "qdec.bin \$$argv" >> $(DESTDIR)$(bindir)/qdec
	chmod a+x $(DESTDIR)$(bindir)/qdec
	if [[ -f /usr/pubsw/bin/chrpath ]]; then \
	/usr/pubsw/bin/chrpath -d $(DESTDIR)$(bindir)/qdec.bin; fi

install-data-local:
	test -z $(DESTDIR)$(foodir) || $(mkdir_p) $(DESTDIR)$(foodir)
		$(INSTALL_DATA) $(foo_datafiles) $(DESTDIR)$(foodir)

uninstall-hook:
	rm -f $(DESTDIR)$(bindir)/qdec.bin

clean-local:
	-rm *Tcl.cxx QdecLibInit.cxx

else
bin_PROGRAMS = 

all:
	@echo "Building KWWidgets apps not enabled."
	@echo "Run configure with either --enable-kww-apps or --with-kwwidgets-dir=DIR"

install-exec-hook:

# these are created on-th-fly (see above), but make dist complains if they
# don't exist, so create some empty files for packaging purposes:
dist-exec-hook:
	@touch vtkKWQdecAppTcl.cxx
	@touch vtkKWQdecWindowTcl.cxx
	@touch vtkKWQdecViewTcl.cxx
	@touch QdecLibInit.cxx
	@sleep 2
	@touch vtkKWQdecApp.h
	@touch vtkKWQdecWindow.h
	@touch vtkKWQdecView.h
	@touch QdecLibInit.data

vtkKWQdecAppTcl.cxx:
vtkKWQdecWindowTcl.cxx:
vtkKWQdecViewTcl.cxx:
QdecLibInit.cxx:
endif

# Copy over precompiled 32-bit versions until we get 64-bit builds working.
if INSTALL_MAC32_GUIS

qdecdir=$(prefix)/lib/resource
qdec_DATA=QdecIcons.txt

install-exec-hook:
	cp $(top_builddir)/distribution/Mac_32b/qdec_Mac32.bin $(DESTDIR)$(bindir)/qdec.bin
	echo "#!/bin/tcsh -ef" > $(DESTDIR)$(bindir)/qdec
	echo "source \$$FREESURFER_HOME/bin/tcl_setup" \
		>> $(DESTDIR)$(bindir)/qdec
	echo "source \$$FREESURFER_HOME/bin/vtk_setup" \
		>> $(DESTDIR)$(bindir)/qdec
	echo "source \$$FREESURFER_HOME/bin/kww_setup" \
		>> $(DESTDIR)$(bindir)/qdec
	echo "if (-e /bin/pwd) then" >> $(DESTDIR)$(bindir)/qdec
	echo "  setenv HOME \`cd \$$HOME; /bin/pwd\`" \
		>> $(DESTDIR)$(bindir)/qdec
	echo "  setenv SUBJECTS_DIR \`cd \$$SUBJECTS_DIR; /bin/pwd\`" \
		>> $(DESTDIR)$(bindir)/qdec
	echo "endif" >> $(DESTDIR)$(bindir)/qdec
	echo "qdec.bin \$$argv" >> $(DESTDIR)$(bindir)/qdec
	chmod a+x $(DESTDIR)$(bindir)/qdec
	if [[ -f /usr/pubsw/bin/chrpath ]]; then \
	/usr/pubsw/bin/chrpath -d $(DESTDIR)$(bindir)/qdec.bin; fi

install-data-local:
	test -z $(DESTDIR)$(foodir) || $(mkdir_p) $(DESTDIR)$(foodir)
		$(INSTALL_DATA) $(foo_datafiles) $(DESTDIR)$(foodir)

uninstall-hook:
	rm -f $(DESTDIR)$(bindir)/qdec.bin
endif

EXTRA_DIST=QdecIcons.txt $(foo_DATA) \
qdec2.table.dat $(dist_bin_SCRIPTS) \
QdecEvents.h \
QdecLibInit.cxx \
QdecLibInit.data \
QdecMain.cxx \
QdecUtilities.cxx \
QdecUtilities.h \
QdecVertexAnnotationLookup.h \
vtkKWQdecApp.cxx \
vtkKWQdecApp.h \
vtkKWQdecAppTcl.cxx \
vtkKWQdecView.cxx \
vtkKWQdecView.h \
vtkKWQdecViewTcl.cxx \
vtkKWQdecWindow.cxx \
vtkKWQdecWindow.h \
vtkKWQdecWindowTcl.cxx \
development_notes.txt \
qdec.ppt

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
