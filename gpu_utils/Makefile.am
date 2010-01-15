##
## Makefile.am 
##

bin_PROGRAMS=cudadetect

cudadetect_SOURCES=cudadetect.cpp
cudadetect_CXXFLAGS=$(CUDA_CFLAGS)
cudadetect_LDFLAGS=-rdynamic -ldl
cudadetect_LDADD=$(CUDA_LIBS)

EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
