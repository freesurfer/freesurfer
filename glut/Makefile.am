#
# Makefile.am
#

AM_CPPFLAGS=$(GL_CFLAGS) -fno-strict-aliasing

if HAVE_OPENGL_LIBS
noinst_LIBRARIES=libglut.a
endif

libglut_a_SOURCES=\
	glut.h \
	glut_8x13.c \
	glut_9x15.c \
	glut_bitmap.c \
	glutbitmap.h \
	glutint.h \
	glut_bwidth.c \
	glut_cindex.c \
	glut_cmap.c \
	glut_cursor.c \
	glut_dials.c \
	glut_dstr.c \
	glut_event.c \
	glut_ext.c \
	glut_fullscrn.c \
	glut_gamemode.c \
	glut_get.c \
	glut_glxext.c \
	glut_hel10.c \
	glut_hel12.c \
	glut_hel18.c \
	glut_init.c \
	glut_input.c \
	glut_joy.c \
	glut_key.c \
	glut_keyctrl.c \
	glut_keyup.c \
	glut_menu.c \
	glut_menu2.c \
	glut_mesa.c \
	glut_modifier.c \
	glut_mroman.c \
	glut_overlay.c \
	glut_roman.c \
	glut_shapes.c \
	glut_space.c \
	glut_stroke.c \
	glutstroke.h \
	glut_swap.c \
	glut_swidth.c \
	glut_tablet.c \
	glut_teapot.c \
	glut_tr10.c \
	glut_tr24.c \
	glut_util.c \
	glut_vidresize.c \
	glut_warp.c \
	glut_win.c \
	glut_winmisc.c \
	layerutil.c \
	layerutil.h

if HAVE_OPENGL_LIBS
check_PROGRAMS=test_glut
endif

test_glut_SOURCES=test_glut.c
test_glut_LDADD=$(addprefix $(top_builddir)/, $(LIBS_MGH)) $(LIBS_VXL) -lz \
	$(LIB_JPEG) $(LIB_TIFF)	$(X_LIBS) $(LIBS_OPENGL) $(top_builddir)/glut/libglut.a

# test_glut app cannot run as auto-test. it requires manual intervention.
# to do that, just uncomment the next line and run 'make check'
# a window should appear where random colored dots appear when moused-over
#TESTS=test_glut

foodir=$(prefix)/docs
foo_DATA=license.glut.txt

EXTRA_DIST=\
	VERSION \
	license.glut.txt

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
