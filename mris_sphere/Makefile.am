## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

bin_PROGRAMS = mris_sphere mris_remove_negative_vertices
mris_sphere_SOURCES=mris_sphere.c
mris_sphere_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_sphere_LDFLAGS=$(OS_LDFLAGS)

mris_remove_negative_vertices_SOURCES=mris_remove_negative_vertices.c
mris_remove_negative_vertices_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_remove_negative_vertices_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
