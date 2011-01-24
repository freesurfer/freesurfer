##
## Makefile.am 
##

bin_PROGRAMS=cudadetect

cudadetect_SOURCES=cudadetect.cpp
cudadetect_LDFLAGS=-rdynamic -ldl

dist_bin_SCRIPTS=cudalib_setup

# BUILDCUDA is defined if configure.in finds CUDA
if BUILDCUDA
cudadetect_CXXFLAGS=$(CUDA_CFLAGS) -DFS_CUDA
cudadetect_LDADD=$(CUDA_LIBS)
endif

# put a wrapper around the bin, to set LD_LIBRARY_PATH to /usr/local/cuda
install-exec-hook:cudadetect
	cp $(top_builddir)/gpu_utils/cudadetect $(DESTDIR)$(bindir)/cudadetect.bin
	echo "#!/bin/tcsh -ef" > $(DESTDIR)$(bindir)/cudadetect
	echo "source \$$FREESURFER_HOME/bin/cudalib_setup" \
		>> $(DESTDIR)$(bindir)/cudadetect
	echo "cudadetect.bin \$$argv" >> $(DESTDIR)$(bindir)/cudadetect
	chmod a+x $(DESTDIR)$(bindir)/cudadetect
	if [[ -f /usr/pubsw/bin/chrpath ]]; then \
	/usr/pubsw/bin/chrpath -d $(DESTDIR)$(bindir)/cudadetect.bin; fi

uninstall-hook:
	rm -f $(DESTDIR)$(bindir)/cudadetect.bin

EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
