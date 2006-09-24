## 
# Makefile.am
##

AM_CPPFLAGS=$(GLUT_CFLAGS) \
	$(TCL_CFLAGS) \
	$(X_CFLAGS) \
	-I$(top_srcdir)/include \
	-I$(top_srcdir)/dicom \
	-I$(top_srcdir)/unix

AM_LDFLAGS=$(LDFLAGS_TCL)

if ENABLE_TCLTK_APPS
bin_PROGRAMS=tksurfer
else
bin_PROGRAMS=
endif
tksurfer_SOURCES=tksurfer.c $(top_srcdir)/fsgdf/fsgdf_wrap.c
tksurfer_LDADD= 	$(top_srcdir)/fsgdf/libfsgdf.a \
  $(addprefix $(top_srcdir)/,$(LIBS_MGH)) \
	$(LIBS_TCL_OPENGL) \
	$(X_LIBS) $(GLUT_LIBS) -lGL
tksurfer_LDFLAGS=

install-exec-hook:
	mv $(bindir)/tksurfer $(bindir)/tksurfer.bin
	cp $(top_srcdir)/scripts/tcltktixblt_wrapper $(bindir)/tksurfer
	echo "tksurfer.bin \$$argv" >> $(bindir)/tksurfer
	chmod a+x $(bindir)/tksurfer

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
