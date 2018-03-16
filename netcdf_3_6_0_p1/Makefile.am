#
# Makefile.am for libnetcdf_3_6_0_p1.a
#

SUBDIRS= 

AM_CPPFLAGS=\
	-I$(top_srcdir)/include -DVERSION=`cat VERSION`

AM_LDFLAGS=

LIBS=

noinst_LIBRARIES=libnetcdf_3_6_0_p1.a

libnetcdf_3_6_0_p1_main_SOURCES=\
attr.c \
dim.c \
error.c \
libvers.c \
nc.c \
ncio.c \
ncx.c \
putget.c \
string.c \
v1hpg.c \
v2i.c \
var.c

libnetcdf_3_6_0_p1_a_SOURCES=\
	$(libnetcdf_3_6_0_p1_main_SOURCES)

G++FLAGS = -Wall -Wextra -Weffc++

install-exec-hook:

uninstall-hook:

EXTRA_DIST=netcdf.h nc.h ncx.h fbits.h rnd.h ncconfig.h ncio.h posixio.c onstack.h

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
