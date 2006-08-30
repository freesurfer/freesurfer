#
# Makefile.am for cephes directory
#

AM_CPPFLAGS=-I$(top_srcdir)/include \
	-DTEEM_DIO=$(TEEM_DIO) -DTEEM_32BIT=$(TEEM_32BIT) \
	-DTEEM_QNANHIBIT=$(TEEM_QNANHIBIT) -DTEEM_ENDIAN=$(TEEM_ENDIAN)

if USE_SC_GSL_REPLACEMENT
noinst_LIBRARIES= libcephes.a
else
noinst_LIBRARIES=
endif
libcephes_a_SOURCES = \
bdtr.c \
chdtr.c \
const.c \
fdtr.c \
gamma.c \
igam.c \
igami.c \
incbet.c \
incbi.c \
mtherr.c \
ndtr.c \
ndtri.c \
polevl.c \
stdtr.c

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
