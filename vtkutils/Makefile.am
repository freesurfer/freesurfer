#
# Makefile.am for libvtkutils.a
#

AM_CPPFLAGS= -Wno-deprecated \
	$(X_CFLAGS) \
	$(TCL_CFLAGS) \
	$(VXL_CFLAGS) \
	-I$(top_srcdir)/include \
	$(VTK_CXXFLAGS) \
	$(KWWIDGETS_CXXFLAGS)

AM_LDFLAGS=

if ENABLE_KWWIDGETS_APPS
noinst_LIBRARIES= libvtkutils.a
else
noinst_LIBRARIES=
endif

libvtkutils_a_SOURCES= \
	vtkKWProgressDialog.cxx vtkKWProgressDialog.h \
	vtkKWRGBATransferFunctionEditor.cxx vtkKWRGBATransferFunctionEditor.h \
	vtkKWRGBATransferFunctionEditorTcl.cxx \
	vtkFreesurferLookupTable.cxx vtkFreesurferLookupTable.h \
	vtkRGBATransferFunction.cxx vtkRGBATransferFunction.h \
	vtkFDTensorGlyph.cxx vtkFDTensorGlyph.h

am_libvtkutils_a_OBJECTS=$(libvtkutils_a_SOURCES:%.cxx=%.o)

# We use vtkWrapTcl to create the Tcl wrapped code for this object.
vtkKWRGBATransferFunctionEditorTcl.cxx: vtkKWRGBATransferFunctionEditor.h
	${VTK_DIR}/bin/vtkWrapTcl vtkKWRGBATransferFunctionEditor.h \
	${VTK_LIB_VTK}/hints 1 vtkKWRGBATransferFunctionEditorTcl.cxx

# .cxx.o:
# 	@echo "compiler: $(COMPILE) end_compiler"
# 	source='$<' object='$@' libtool=no \
# 	depfile='$(DEPDIR)/$*.Po' tmpdepfile='$(DEPDIR)/$*.TPo' \
# 	$(CCDEPMODE) $(depcomp) \
# 	$(COMPILE) -c `test -f '$<' || echo '$(srcdir)/'`$<

EXTRA_DIST=$(wildcard *.h)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
