## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

bin_PROGRAMS = mris_fix_topology mris_remove_intersection

mris_fix_topology_SOURCES=mris_fix_topology.c
mris_fix_topology_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_fix_topology_LDFLAGS=$(OS_LDFLAGS)

mris_remove_intersection_SOURCES=mris_remove_intersection.c
mris_remove_intersection_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_remove_intersection_LDFLAGS=$(OS_LDFLAGS)

foodir=$(prefix)/docs/xml
foo_DATA=mris_fix_topology.help.xml mris_remove_intersection.help.xml

## ----
## CUDA
## ----

# BUILDCUDA is defined if configure.in finds CUDA
if BUILDCUDA

NVCCFLAGS += -g -arch sm_13 --ptxas-options=-v

# rules for building cuda files
.cu.o:
	$(NVCC) -o $@ -c $< $(NVCCFLAGS) $(CUDA_CFLAGS) -I$(top_srcdir)/include

bin_PROGRAMS += mris_fix_topology_cuda
mris_fix_topology_cuda_SOURCES = mris_fix_topology.c
mris_fix_topology_cuda_CFLAGS = $(CUDA_CFLAGS) $(AM_CFLAGS) -DFS_CUDA
mris_fix_topology_cuda_CXXFLAGS = $(CUDA_CFLAGS) $(AM_CFLAGS) -DFS_CUDA
mris_fix_topology_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS)
mris_fix_topology_cuda_LDFLAGS = $(OS_LDFLAGS) 
mris_fix_topology_cuda_LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(mris_fix_topology_cuda_CFLAGS) \
        $(CFLAGS) $(mris_fix_topology_cuda_LDFLAGS) $(LDFLAGS) -o $@

endif

TESTS=$(top_builddir)/scripts/help_xml_validate

EXTRA_DIST=$(foo_DATA) $(BUILT_SOURCES)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra

BUILT_SOURCES=mris_fix_topology.help.xml.h mris_remove_intersection.help.xml.h
