##
## Makefile.am 
##

bin_PROGRAMS=cudadetect

cudadetect_SOURCES=cudadetect.cpp
cudadetect_LDFLAGS=-rdynamic -ldl

# BUILDCUDA is defined if configure.in finds CUDA
if BUILDCUDA
cudadetect_CXXFLAGS=$(CUDA_CFLAGS) -DFS_CUDA
cudadetect_LDADD=$(CUDA_LIBS)
endif

EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
