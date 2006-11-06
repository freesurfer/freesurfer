## 
# Makefile.am
##

AM_CFLAGS=$(GLUT_CFLAGS) \
	-I$(top_srcdir)/include \
	-I$(top_srcdir)/include/dicom \
	$(X_CFLAGS)

if HAVE_GLUT_LIBS
bin_PROGRAMS = test_window_env
test_window_env_SOURCES=test_window_env.c
test_window_env_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH)) \
	$(GLUT_LIBS) \
	$(X_LIBS) -lGL  
test_window_env_LDFLAGS=
else
bin_PROGRAMS=
endif

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
