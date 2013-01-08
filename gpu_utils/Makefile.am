##
## Makefile.am 
##

# BUILDCUDA is defined if configure.in finds CUDA
if BUILDCUDA

bin_PROGRAMS=cudadetect

cudadetect_SOURCES=cudadetect.cpp


if HAVE_MAC_OSX
cudadetect_LDFLAGS=-rdynamic -ldl -F/Library/Frameworks -framework CUDA
else
cudadetect_LDFLAGS=-rdynamic -ldl
endif

cudadetect_CXXFLAGS=$(CUDA_CFLAGS) -DFS_CUDA
cudadetect_LDADD=$(CUDA_LIBS)

# put a wrapper around the bin, to set LD_LIBRARY_PATH to /usr/local/cuda
install-exec-hook:cudadetect
	cp $(top_builddir)/gpu_utils/cudadetect $(DESTDIR)$(bindir)/cudadetect.bin
	echo "#!/bin/tcsh -ef" > $(DESTDIR)$(bindir)/cudadetect
	echo "source \$$FREESURFER_HOME/bin/cuda_setup" \
		>> $(DESTDIR)$(bindir)/cudadetect
	echo "cudadetect.bin \$$argv" >> $(DESTDIR)$(bindir)/cudadetect
	chmod a+x $(DESTDIR)$(bindir)/cudadetect
	if [[ -f /usr/pubsw/bin/chrpath ]]; then \
	/usr/pubsw/bin/chrpath -d $(DESTDIR)$(bindir)/cudadetect.bin; fi

uninstall-hook:
	rm -f $(DESTDIR)$(bindir)/cudadetect.bin

endif


EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
