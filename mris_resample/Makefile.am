## 
## Makefile.am 
##

if HAVE_ANN_LIBS
if HAVE_BOOST_LIBS

AM_CXXFLAGS=-I$(top_srcdir)/include $(ANN_CFLAGS) $(BOOST_CFLAGS)

bin_PROGRAMS = mris_resample
mris_resample_SOURCES=mris_resample.cpp
mris_resample_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH)) \
	$(BOOST_PO_LIB) $(ANN_LIBS) -lANN
mris_resample_LDFLAGS=$(OS_LDFLAGS)

endif
endif

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
