## 
## Makefile.am 
##

AM_CPPFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

check_PROGRAMS = test_mriSoapBubbleFloat

TESTS=test_mriSoapBubbleFloat_withData

test_mriSoapBubbleFloat_SOURCES=test_mriSoapBubbleFloat.cpp
test_mriSoapBubbleFloat_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
test_mriSoapBubbleFloat_LDFLAGS= $(OS_LDFLAGS)

EXTRA_DIST=test_mriSoapBubbleFloat_withData testdata.tar.gz

EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra

clean-local:
	rm -f *.o
