#
# Makefile.am
#

AM_CPPFLAGS=-fno-strict-aliasing -DHAVE_CONFIG_H -D_REENTRANT

AM_LDFLAGS=

LIBS=

noinst_LIBRARIES=libxml2.a

libxml2_a_SOURCES=\
	c14n.c \
	catalog.c \
	chvalid.c \
	debugXML.c \
	dict.c \
	DOCBparser.c \
	encoding.c \
	entities.c \
	error.c \
	globals.c \
	hash.c \
	HTMLparser.c \
	HTMLtree.c \
	list.c \
	nanoftp.c \
	nanohttp.c \
	parser.c \
	parserInternals.c \
	pattern.c \
	relaxng.c \
	SAX2.c \
	SAX.c \
	schematron.c \
	threads.c \
	tree.c \
	uri.c \
	valid.c \
	xinclude.c \
	xlink.c \
	xmlIO.c \
	xmlmemory.c \
	xmlmodule.c \
	xmlreader.c \
	xmlregexp.c \
	xmlsave.c \
	xmlschemas.c \
	xmlschemastypes.c \
	xmlstring.c \
	xmlunicode.c \
	xmlwriter.c \
	xpath.c \
	xpointer.c \
	libxml/c14n.h \
	libxml/catalog.h \
	libxml/chvalid.h \
	libxml/debugXML.h \
	libxml/dict.h \
	libxml/DOCBparser.h \
	libxml/encoding.h \
	libxml/entities.h \
	libxml/globals.h \
	libxml/hash.h \
	libxml/HTMLparser.h \
	libxml/HTMLtree.h \
	libxml/list.h \
	libxml/nanoftp.h \
	libxml/nanohttp.h \
	libxml/parser.h \
	libxml/parserInternals.h \
	libxml/pattern.h \
	libxml/relaxng.h \
	libxml/SAX2.h \
	libxml/SAX.h \
	libxml/schemasInternals.h \
	libxml/schematron.h \
	libxml/threads.h \
	libxml/tree.h \
	libxml/uri.h \
	libxml/valid.h \
	libxml/xinclude.h \
	libxml/xlink.h \
	libxml/xmlautomata.h \
	libxml/xmlerror.h \
	libxml/xmlexports.h \
	libxml/xmlIO.h \
	libxml/xmlmemory.h \
	libxml/xmlmodule.h \
	libxml/xmlreader.h \
	libxml/xmlregexp.h \
	libxml/xmlsave.h \
	libxml/xmlschemas.h \
	libxml/xmlschemastypes.h \
	libxml/xmlstring.h \
	libxml/xmlunicode.h \
	libxml/xmlwriter.h \
	libxml/xpath.h \
	libxml/xpathInternals.h \
	libxml/xpointer.h

check_PROGRAMS=runtest testapi testchar testdict testrecurse

runtest_SOURCES=runtest.c
runtest_LDADD=$(top_builddir)/xml2/libxml2.a -ldl -lpthread -lz -lm

testapi_SOURCES=testapi.c legacy.c
testapi_LDADD=$(top_builddir)/xml2/libxml2.a -ldl -lpthread -lz -lm

testchar_SOURCES=testchar.c
testchar_LDADD=$(top_builddir)/xml2/libxml2.a -lpthread -lz -lm

testdict_SOURCES=testdict.c
testdict_LDADD=$(top_builddir)/xml2/libxml2.a -lz

testrecurse_SOURCES=testrecurse.c
testrecurse_LDADD=$(top_builddir)/xml2/libxml2.a -lz

TESTS=runtest testapi testchar testdict testrecurse

foodir=$(prefix)/docs
foo_DATA=license.xml2

EXTRA_DIST=trionan.c \
	README \
	license.libxml2

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
