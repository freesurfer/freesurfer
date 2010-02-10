## 
# Makefile.am
##

AM_CFLAGS=-I$(top_srcdir)/include $(GL_CFLAGS)

if HAVE_OPENGL_LIBS
check_PROGRAMS = test_window_env
test_window_env_SOURCES=test_window_env.c
test_window_env_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH)) \
	$(LIBS_OPENGL) $(X_LIBS) -lGL
test_window_env_LDFLAGS=
if HAVE_MAC_OSX
TESTS=
else
TESTS=test_window_env
endif
else
bin_PROGRAMS=
endif

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
