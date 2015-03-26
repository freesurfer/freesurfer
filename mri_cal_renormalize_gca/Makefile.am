## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include

AM_CXXFLAGS=-I$(top_srcdir)/include \
             $(BOOST_CFLAGS)

AM_LDFLAGS=

bin_PROGRAMS = mri_cal_renormalize_gca

mri_cal_renormalize_gca_SOURCES=mri_cal_renormalize_gca.c
mri_cal_renormalize_gca_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH)) 
mri_cal_renormalize_gca_LDFLAGS=$(OS_LDFLAGS)

#foodir=$(prefix)/docs/xml
#foo_DATA=mri_cal_renormalize_gca.help.xml mri_cal_renormalize_gca_cuda.help.xml
#foo2dir=$(prefix)/docs/html
#foo2_DATA=mri_cal_renormalize_gca.help.xml.html mri_cal_renormalize_gca_cuda.help.xml.html

## ----
## CUDA
## ----


TESTS=$(top_builddir)/scripts/help_xml_validate

clean-local:
	rm -f $(BUILT_SOURCES)

EXTRA_DIST=$(foo_DATA) $(BUILT_SOURCES)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra

BUILT_SOURCES=mri_cal_renormalize_gca.help.xml.h  \
	mri_cal_renormalize_gca.help.xml.html 
