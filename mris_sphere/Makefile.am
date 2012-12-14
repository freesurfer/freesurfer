## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

bin_PROGRAMS = mris_sphere mris_remove_negative_vertices

mris_sphere_SOURCES=mris_sphere.c mris_sphere.help.xml.h
mris_sphere_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_sphere_LDFLAGS=$(OS_LDFLAGS)

foodir=$(prefix)/docs/xml
foo_DATA=mris_sphere.help.xml

mris_remove_negative_vertices_SOURCES=mris_remove_negative_vertices.c
mris_remove_negative_vertices_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_remove_negative_vertices_LDFLAGS=$(OS_LDFLAGS)

## ----
## CUDA
## ----

# BUILDCUDA is defined if configure.in finds CUDA
if BUILDCUDA

NVCCFLAGS += -g -arch sm_13 --ptxas-options=-v

# rules for building cuda files
.cu.o:
	$(NVCC) -o $@ -c $< $(NVCCFLAGS) $(CUDA_CFLAGS) -I$(top_srcdir)/include

bin_PROGRAMS += mris_sphere_cuda
mris_sphere_cuda_SOURCES = mris_sphere.c
mris_sphere_cuda_CFLAGS = $(CUDA_CFLAGS) $(AM_CFLAGS) -DFS_CUDA
mris_sphere_cuda_CXXFLAGS = $(CUDA_CFLAGS) $(AM_CFLAGS) -DFS_CUDA
mris_sphere_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS)
mris_sphere_cuda_LDFLAGS = $(OS_LDFLAGS) 
mris_sphere_cuda_LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(mris_sphere_cuda_CFLAGS) \
        $(CFLAGS) $(mris_sphere_cuda_LDFLAGS) $(LDFLAGS) -o $@

endif

TESTS=$(top_builddir)/scripts/help_xml_validate

EXTRA_DIST=$(foo_DATA) $(BUILT_SOURCES) make_mris_sphere_movie.tcl

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra

BUILT_SOURCES=mris_sphere.help.xml.h

