#
# Makefile.am
#

AM_CPPFLAGS=

AM_LDFLAGS=

LIBS=

noinst_LIBRARIES=libjpeg.a

libjpeg_a_SOURCES=\
	jcapimin.c \
	jcapistd.c \
	jctrans.c \
	jcparam.c \
	jdatadst.c \
	jcinit.c \
	jcmaster.c \
	jcmarker.c \
	jcmainct.c \
	jcprepct.c \
	jccoefct.c \
	jccolor.c \
	jcsample.c \
	jchuff.c \
	jcphuff.c \
	jcdctmgr.c \
	jfdctfst.c \
	jfdctflt.c \
	jfdctint.c \
	jdapimin.c \
	jdapistd.c \
	jdtrans.c \
	jdatasrc.c \
	jdmaster.c \
	jdinput.c \
	jdmarker.c \
	jdhuff.c \
	jdphuff.c \
	jdmainct.c \
	jdcoefct.c \
	jdpostct.c \
	jddctmgr.c \
	jidctfst.c \
	jidctflt.c \
	jidctint.c \
	jidctred.c \
	jdsample.c \
	jdcolor.c \
	jquant1.c \
	jquant2.c \
	jdmerge.c \
	jcomapi.c \
	jutils.c \
	jerror.c \
	jmemmgr.c \
	jmemnobs.c

cjpeg_SOURCES=\
	cjpeg.c \
	rdppm.c \
	rdgif.c \
	rdtarga.c \
	rdrle.c \
	rdbmp.c \
	rdswitch.c \
	cdjpeg.c

cjpeg_LDADD=libjpeg.a

djpeg_SOURCES=\
	djpeg.c \
	wrppm.c \
	wrgif.c \
	wrtarga.c \
	wrrle.c \
	wrbmp.c \
	rdcolmap.c \
	cdjpeg.c

djpeg_LDADD=libjpeg.a

jpegtran_SOURCES=\
	jpegtran.c \
	rdswitch.c \
	cdjpeg.c \
	transupp.c

jpegtran_LDADD=libjpeg.a

rdjpgcom_SOURCES=rdjpgcom.c

wrjpgcom_SOURCES=wrjpgcom.c

# trick to get test data into the build directory
foo_SOURCES=
foo:
	if ! test -f "$(top_builddir)/jpeg/testorig.jpg"; then \
		cp -v $(top_srcdir)/jpeg/test* \
					$(top_builddir)/jpeg/; fi;

check_PROGRAMS=foo cjpeg djpeg jpegtran rdjpgcom wrjpgcom

TESTS=test_jpeg

clean-local:
	rm -f testout*

foodir=$(prefix)/docs
foo_DATA=license.jpeg.txt

EXTRA_DIST=\
	cderror.h  \
	cdjpeg.h  \
	jchuff.h  \
	jconfig.h  \
	jdct.h  \
	jdhuff.h  \
	jerror.h  \
	jinclude.h  \
	jmemsys.h  \
	jmorecfg.h  \
	jpegint.h  \
	jpeglib.h  \
	jversion.h  \
	transupp.h \
	coderules.doc \
	filelist.doc  \
	install.doc  \
	jconfig.doc  \
	libjpeg.doc  \
	structure.doc  \
	usage.doc  \
	wizard.doc \
	README \
	license.jpeg.txt \
	test_jpeg \
	testorig.jpg \
	testimg.ppm \
	testprog.jpg \
	testimg.bmp \
	testimg.jpg \
	testimgp.jpg

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
