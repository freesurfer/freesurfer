##
## Makefile.am 
##
#######################################################################
# When you create your own program directory, copy the 'from' to the 
# 'to' lines to your Makefile.am and change dummy to your program name.
#  
# 0. mkdir  dev/yourdir and put .c file in there
# 1. modify dev/configure.in to have yourdir/Makefile
# 2. modify dev/Makefile.am to have yourdir as one of SUBDIRS
# 3. create dev/yourdir/Makefile.am  (follow dummy/Makefile.am)
#      i.e. replace dummy with your program name.
#
# 4. automake -v at the top directory, i.e. dev/ (which will create 
#    Makefile.in in yourdir directory)
# 5. autoconf (which will create configure from configure.in)
# 6. configure .... and make to make sure that your new thing is working
# 7. try make to make sure that yourthing is compiled.
#
# After confirming everything working:
# 8. checkin yourdir/Makefile.am, yourdir/yourprog.c or cpp.  
#    (Don't checkin Makefile or Makefile.in.)
# 9. checkin dev/configure.in, dev/configure, dev/Makefile.am
#
# If you have problems, please let us know (kteich@nmr.mgh.harvard.edu or
# nicks@nmr.mgh.harvard.edu).
#########################################################################
# from ------------------------------------------------
AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mris_seg2annot
mris_seg2annot_SOURCES=mris_seg2annot.c
mris_seg2annot_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_seg2annot_LDFLAGS=$(OS_LDFLAGS)
# to --------------------------------------------------
#
# Where $(OS_LDFLAGS) is to compile programs statically for Linux.
#
# You can have multiple bin_PROGRAMS in the same directory.
#

EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
