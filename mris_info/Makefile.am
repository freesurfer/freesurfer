## 
## Makefile.am 
##

if USE_LOCAL_EXPAT
EXPAT_THINGS=-I$(top_srcdir)/expat
else
EXPAT_THINGS=
endif

AM_CXXFLAGS=-I$(top_srcdir)/include $(EXPAT_THINGS)
AM_LDFLAGS=

bin_PROGRAMS = mris_info
mris_info_SOURCES=mris_info.cpp
mris_info_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_info_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
