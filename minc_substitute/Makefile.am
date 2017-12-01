#
# Makefile.am for libminc_substitute.a
#

SUBDIRS= 

AM_CPPFLAGS=\
	-I$(top_srcdir)/include -DHAVE_ZLIB

AM_LDFLAGS=

LIBS=

noinst_LIBRARIES=libminc_substitute.a

libminc_substitute_main_SOURCES=\
	minc_convenience.c \
	minc_files.c \
	minc_gaussian.c \
	minc_image_conversion.c \
	minc_input_free.c \
	minc_input_mnc.c \
	minc_inverse.c \
	minc_multidim_arrays.c \
	minc_netcdf_convenience.c \
	minc_output_mnc.c \
	minc_output_volume.c \
	minc_progress.c \
	minc_time.c \
	minc_value_conversion.c \
	minc_volume_io.c

libminc_substitute_a_SOURCES=\
	$(libminc_substitute_main_SOURCES)

G++FLAGS = -Wall -Wextra -Weffc++

install-exec-hook:

uninstall-hook:

EXTRA_DIST=

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
