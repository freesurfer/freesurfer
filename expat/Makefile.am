#
# Makefile.am
#

AM_CPPFLAGS=-fexceptions

AM_CXXFLAGS=-fexceptions

AM_LDFLAGS=

LIBS=

noinst_LIBRARIES=libexpat.a

libexpat_a_SOURCES=\
	xmlparse.c ascii.h internal.h \
	xmltok.c xmltok.h xmltok_impl.h \
	nametab.h asciitab.h utf8tab.h iasciitab.h latin1tab.h \
	xmlrole.c xmlrole.h \
	xmltchar.h \
	xmlfile.c xmlfile.h \
	codepage.c codepage.h \
	unixfilemap.c filemap.h \
	expat.h expat_external.h

check_PROGRAMS=runtests runtestspp

runtests_SOURCES=runtests.c chardata.c chardata.h minicheck.c minicheck.h
runtests_LDADD=$(top_builddir)/expat/libexpat.a

runtestspp_SOURCES=runtestspp.cpp chardata.c chardata.h minicheck.c minicheck.h
runtestspp_LDADD=$(top_builddir)/expat/libexpat.a

TESTS=runtests runtestspp

foodir=$(prefix)/docs
foo_DATA=license.expat.txt

EXTRA_DIST=\
	xmltok_impl.c \
	xmltok_ns.c \
	README \
	license.expat.txt

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
