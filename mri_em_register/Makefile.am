## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include

AM_CXXFLAGS=-I$(top_srcdir)/include\
             $(BOOST_CFLAGS)\
             -Wall -Wextra

AM_LDFLAGS=

bin_PROGRAMS = mri_em_register
mri_em_register_SOURCES=mri_em_register.c\
                        findtranslation.cpp findtranslation.h\
                        emregisterutils.cpp emregisterutils.h
mri_em_register_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_em_register_LDFLAGS=$(OS_LDFLAGS)


foodir=$(prefix)/docs/xml
foo_DATA=mri_em_register.help.xml mri_em_register_cuda.help.xml


## ----
## CUDA
## ----

# BUILDCUDA is defined if configure.in finds CUDA
if BUILDCUDA
# rules for building cuda files
.cu.o:
	$(NVCC) -o $@ -c $< $(NVCCFLAGS) $(CUDA_CFLAGS) $(AM_CFLAGS) $(MNI_CFLAGS)

bin_PROGRAMS += mri_em_register_cuda

mri_em_register_cuda_SOURCES = mri_em_register.c \
  em_register_cuda.cu em_register_cuda.h\
  findtranslation.cpp findtranslation.h\
  emregisterutils.cpp emregisterutils.h
mri_em_register_cuda_CFLAGS = $(CUDA_CFLAGS) $(AM_CFLAGS) -DFS_CUDA
mri_em_register_cuda_CXXFLAGS = $(CUDA_CFLAGS) $(AM_CFLAGS) -DFS_CUDA
mri_em_register_cuda_LDADD = $(addprefix $(top_builddir)/, $(LIBS_CUDA_MGH)) $(CUDA_LIBS)
mri_em_register_cuda_LDFLAGS = $(OS_LDFLAGS) 
mri_em_register_cuda_LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(mri_em_register_cuda_CFLAGS) \
        $(CFLAGS) $(mri_em_register_cuda_LDFLAGS) $(LDFLAGS) -o $@
endif

TESTS=$(top_builddir)/scripts/help_xml_validate

EXTRA_DIST=$(foo_DATA) $(BUILT_SOURCES)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra

BUILT_SOURCES=mri_em_register.help.xml.h mri_em_register_cuda.help.xml.h
