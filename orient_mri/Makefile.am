AM_CFLAGS=-I$(top_srcdir)/include \
	$(VTK_CFLAGS) $(KWWIDGETS_CFLAGS)
AM_CXXFLAGS=-Wno-deprecated -I$(top_srcdir)/include \
	$(VTK_CXXFLAGS) $(KWWIDGETS_CXXFLAGS)

if ENABLE_KWWIDGETS_APPS
bin_PROGRAMS = orient_mri
orient_mri_SOURCES= \
	OrientMRIMain.cxx OrientMRIMain.h \
	vtkKWOrientMRIApp.cxx vtkKWOrientMRIApp.h vtkKWOrientMRIAppTcl.cxx \
	vtkKWOrientMRIWindow.cxx vtkKWOrientMRIWindow.h vtkKWOrientMRIWindowTcl.cxx \
	vtkKWOrientMRIView.cxx vtkKWOrientMRIView.h vtkKWOrientMRIViewTcl.cxx \
	vtkFSVolumeSource.cxx vtkFSVolumeSource.h \
	OrientMRILibInit.cxx
orient_mri_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH)) \
	$(VTK_LIBS) -lvtkWidgetsTCL -lvtkWidgets -lvtkHybridTCL \
	-lvtkHybrid -lvtkexoIIc -lvtkNetCDF -lvtkVolumeRenderingTCL \
	-lvtkVolumeRendering -lvtkRenderingTCL -lvtkRendering -lvtkftgl \
	-lvtkfreetype  -L/tmp $(KWWIDGETS_LIBS) -L/tmp
orient_mri_LDFLAGS=$(OS_LDFLAGS)

vtkKWOrientMRIAppTcl.cxx: vtkKWOrientMRIApp.h
	${VTK_DIR}/bin/vtkWrapTcl vtkKWOrientMRIApp.h \
	${VTK_DIR}/src/VTK/Wrapping/hints 1 vtkKWOrientMRIAppTcl.cxx

vtkKWOrientMRIWindowTcl.cxx: vtkKWOrientMRIWindow.h
	${VTK_DIR}/bin/vtkWrapTcl vtkKWOrientMRIWindow.h \
	${VTK_DIR}/src/VTK/Wrapping/hints 1 vtkKWOrientMRIWindowTcl.cxx

vtkKWOrientMRIViewTcl.cxx: vtkKWOrientMRIView.h
	${VTK_DIR}/bin/vtkWrapTcl vtkKWOrientMRIView.h \
	${VTK_DIR}/src/VTK/Wrapping/hints 1 vtkKWOrientMRIViewTcl.cxx

OrientMRILibInit.cxx: OrientMRILibInit.data
	${VTK_DIR}/bin/vtkWrapTclInit OrientMRILibInit.data OrientMRILibInit.cxx

else
bin_PROGRAMS = 
endif

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
