#
# Makefile.am
#
# tiffconf.h defines the build options (JPEG_SUPPORT, ZIP_SUPPORT)

if USE_LOCAL_JPEG
JPEGCFLAGS=-I$(top_srcdir)/jpeg
LIBJPEG=$(top_builddir)/jpeg/libjpeg.a
else
JPEGCFLAGS=$(JPEG_CFLAGS)
LIBJPEG=$(LIB_JPEG)
endif

AM_CPPFLAGS=$(JPEGCFLAGS)

AM_LDFLAGS=

LIBS=

noinst_LIBRARIES=libtiff.a

libtiff_a_SOURCES=\
	uvcode.h \
	port.h \
	tiffcomp.h \
	tiff.h \
	tiffiop.h \
	t4.h \
	tiffconf.h \
	tiffio.h \
	tiffvers.h \
	tif_aux.c \
	tif_close.c \
	tif_codec.c \
	tif_color.c \
	tif_compress.c \
	tif_dir.c \
	tif_dir.h \
	tif_dirinfo.c \
	tif_dirread.c \
	tif_dirwrite.c \
	tif_dumpmode.c \
	tif_error.c \
	tif_extension.c \
	tif_flush.c \
	tif_getimage.c \
	tif_jpeg.c \
	tif_luv.c \
	tif_lzw.c \
	tif_next.c \
	tif_ojpeg.c \
	tif_open.c \
	tif_packbits.c \
	tif_pixarlog.c \
	tif_predict.c \
	tif_predict.h \
	tif_print.c \
	tif_read.c \
	tif_strip.c \
	tif_swab.c \
	tif_thunder.c \
	tif_tile.c \
	tif_unix.c \
	tif_version.c \
	tif_warning.c \
	tif_write.c \
	tif_zip.c

# trick to get test data into build dir
foo_SOURCES=
foo:
	if ! test -f "$(top_builddir)/tiff/icon_folder.tif"; then \
		cp -v $(top_srcdir)/tiff/icon_folder.tif \
					$(top_builddir)/tiff/; fi;

check_PROGRAMS=foo test_tiff

test_tiff_SOURCES=test_tiff.c
test_tiff_LDADD=$(top_builddir)/tiff/libtiff.a \
	$(LIBJPEG) -lz

TESTS=tiff_test

foodir=$(prefix)/docs
foo_DATA=license.tiff.txt

EXTRA_DIST=\
	tiff_test \
	icon_folder.tif \
	VERSION \
	license.tiff.txt

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
