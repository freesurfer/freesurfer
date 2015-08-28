##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mri_correct_segmentations
mri_correct_segmentations_SOURCES=mri_correct_segmentations.c
mri_correct_segmentations_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_correct_segmentations_LDFLAGS=$(OS_LDFLAGS)

#
# Where $(OS_LDFLAGS) is to compile programs statically for Linux.
#
# You can have multiple bin_PROGRAMS in the same directory.
#
########################################################################
# When you want to create your program to test, use the following 
# framework:
# You can have as many programs (separated by space) in check_PROGRAMS.
# Your program can consist of many sources (you add to _SOURCES= line).
# You can do "make myown" to compile. The default CFLAGS="-g -O2", but
# you can do "make CFLAGS=-g myown" to create a debugging version.
#
# The program listed as check_PROGRAMS are not compiled under standard 
# "make". Only "make check" will compile them.
# 
# Please don't remove myown.c from this mri_correct_segmentations directory.
########################################################################
#check_PROGRAMS=myown

#myown_SOURCES=myown.c
#myown_LDADD=$(addprefix $(top_builddir)/, $(LIBS_MGH))

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
