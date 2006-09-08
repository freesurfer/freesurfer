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
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom $(VTK_CXXFLAGS)/vtk-5.0 -Wno-deprecated

bin_PROGRAMS = mris_volmask
mris_volmask_SOURCES=mris_volmask.cpp cmd_line_interface.cpp
mris_volmask_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH)) -lvtkFiltering -lvtkHybrid -lvtkRendering -lvtkGraphics -lvtkImaging -lvtkftgl -lvtkfreetype -lGL -lXt -lSM -lICE -lSM -lICE -lX11 -lXext -lX11 -lXext -lvtkIO -lvtkFiltering -lvtkCommon -lpthread -lm -lvtkDICOMParser -lvtkpng -lvtktiff -lvtkzlib -lvtkjpeg -lvtkexpat -lvtksys -ldl -lvtkMPEG2Encode -lvtkexoIIc -lvtkNetCDF
mris_volmask_LDFLAGS=$(OS_LDFLAGS) $(VTK_LIBS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
