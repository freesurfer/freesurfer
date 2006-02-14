#
# Makefile.am for NrrdIO directory
#

AM_CPPFLAGS=-I$(top_srcdir)/include/NrrdIO \
	-DTEEM_DIO=$(TEEM_DIO) -DTEEM_32BIT=$(TEEM_32BIT) \
	-DTEEM_QNANHIBIT=$(TEEM_QNANHIBIT) -DTEEM_ENDIAN=$(TEEM_ENDIAN)

noinst_LIBRARIES= libnrrdio.a
libnrrdio_a_SOURCES = 754.c mop.c array.c parseAir.c dio.c sane.c endianAir.c string.c enum.c miscAir.c biffbiff.c accessors.c defaultsNrrd.c enumsNrrd.c arraysNrrd.c methodsNrrd.c reorder.c axis.c simple.c comment.c keyvalue.c endianNrrd.c parseNrrd.c gzio.c read.c write.c format.c formatNRRD.c encoding.c encodingRaw.c encodingAscii.c encodingHex.c encodingGzip.c subset.c encodingBzip2.c formatEPS.c formatPNG.c formatPNM.c formatText.c formatVTK.c


# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
