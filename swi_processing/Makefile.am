## 
## Makefile.am 
##

AM_CXXFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

bin_PROGRAMS = swi_preprocess swi_process
swi_preprocess_SOURCES=swi_preprocess.cpp cmd_line_interface.cpp
swi_preprocess_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
swi_preprocess_LDFLAGS=$(OS_LDFLAGS)

swi_process_SOURCES=swi_process.cpp cmd_line_interface.cpp
swi_process_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
swi_process_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
