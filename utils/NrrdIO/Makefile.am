#
# Makefile.am for NrrdIO directory
#

AM_CPPFLAGS=-I$(top_srcdir)/include/NrrdIO \
	-DTEEM_DIO=$(TEEM_DIO) -DTEEM_32BIT=$(TEEM_32BIT) \
	-DTEEM_QNANHIBIT=$(TEEM_QNANHIBIT) -DTEEM_ENDIAN=$(TEEM_ENDIAN)

noinst_LIBRARIES= libnrrdio.a

libnrrdio_a_SOURCES = \
	754.c            formatPNG.c \
	accessors.c      formatPNM.c \
	array.c          formatText.c \
	arraysNrrd.c     formatVTK.c \
	axis.c           gzio.c \
	biffbiff.c       keyvalue.c \
	comment.c        methodsNrrd.c \
	defaultsNrrd.c   miscAir.c \
	dio.c            mop.c \
	encodingAscii.c  parseAir.c \
	encodingBzip2.c  parseNrrd.c \
	encoding.c       preamble.c \
	encodingGzip.c   qnanhibit.c \
	encodingHex.c    read.c \
	encodingRaw.c    reorder.c \
	endianAir.c      sampleIO.c \
	endianNrrd.c     sane.c \
	enum.c           simple.c \
	enumsNrrd.c      string.c \
	format.c         subset.c \
	formatEPS.c      write.c \
	formatNRRD.c

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
