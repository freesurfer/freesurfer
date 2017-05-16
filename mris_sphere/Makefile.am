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
foo2dir=$(prefix)/docs/html
foo2_DATA=mris_sphere.help.xml.html

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

# trick to get test data into the build directory
foo_SOURCES=
foo:
	if ! test -f "$(top_builddir)/mris_sphere/testdata.tar.gz"; then \
		cp -v $(top_srcdir)/mris_sphere/testdata.tar.gz \
					$(top_builddir)/mris_sphere/; fi;

check_PROGRAMS=foo

TESTS=test_mris_sphere test_fs_posix $(top_builddir)/scripts/help_xml_validate

clean-local:
	rm -f $(BUILT_SOURCES)

EXTRA_DIST=$(foo_DATA) $(BUILT_SOURCES) \
	make_mris_sphere_movie.tcl test_mris_sphere

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra

BUILT_SOURCES=mris_sphere.help.xml.h mris_sphere.help.xml.html

