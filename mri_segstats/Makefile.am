## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mri_segstats

mri_segstats_SOURCES=mri_segstats.c mri_segstats.help.xml.h
mri_segstats_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_segstats_LDFLAGS=$(OS_LDFLAGS)


## ----
## CUDA
## ----

# BUILDCUDA is defined if configure.in finds CUDA
if BUILDCUDA

# rules for building cuda files
.cu.o:
	$(NVCC) -o $@ -c $< $(NVCCFLAGS) $(CUDA_CFLAGS) $(AM_CFLAGS) $(MNI_CFLAGS)

bin_PROGRAMS += mri_segstats_cuda

mri_segstats_cuda_SOURCES = mri_segstats.c
mri_segstats_cuda_CFLAGS = $(CUDA_CFLAGS) $(AM_CFLAGS) -DFS_CUDA
mri_segstats_cuda_CXXFLAGS = $(CUDA_CFLAGS) $(AM_CFLAGS) -DFS_CUDA
mri_segstats_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS)
mri_segstats_cuda_LDFLAGS = $(OS_LDFLAGS) 
mri_segstats_cuda_LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(mri_segstats_cuda_CFLAGS) \
        $(CFLAGS) $(mri_segstats_cuda_LDFLAGS) $(LDFLAGS) -o $@

endif
# endif BUILDCUDA



foodir=$(prefix)/docs/xml
foo_DATA=mri_segstats.help.xml

# trick to get test data into the build directory
foo_SOURCES=
foo:
	if ! test -f "$(top_builddir)/mri_segstats/testdata.tar.gz"; then \
		cp -v $(top_srcdir)/mri_segstats/testdata.tar.gz \
					$(top_builddir)/mri_segstats/; fi;

check_PROGRAMS=foo

TESTS=test_mri_segstats $(top_builddir)/scripts/help_xml_validate

EXTRA_DIST=test_mri_segstats testdata.tar.gz $(foo_DATA) $(BUILT_SOURCES)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra

BUILT_SOURCES=mri_segstats.help.xml.h

