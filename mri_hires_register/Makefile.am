## 
## Makefile.am 
##
#######################################################################
# when you create your own program directory, copy the from to the to lines
# to your Makefile.am and change mri_hires_register to your program name
#  
# 0. mkdir  dev/yourdir and put .c file in there                   <--- old way
# 1. modify dev/configure.in to have yourdir/Makefile              
# 2. modify dev/Makefile.am to have yourdir as one of SUBDIRS      
# 3. create dev/yourdir/Makefile.am  (follow mri_hires_register/Makefile.am)    
#      i.e. replace mri_hires_register with your program name.
#
# 4. automake -v at the top directory, i.e. dev/ (which will create Makefile.in in yourthing directory)
# 4. autoconf (which will create configure from configure.in)
# 5. configure .... and make to make sure that your new thing is working
# 6. try make to make sure that yourthing is compiled.
#
# After confirming everything working 
# 6. checkin yourdir/Makefile.am, yourdir/Makefile.in, yourdir/yourprog.c or cpp.  (Don't checkin Makefile.)
# 7. checkin dev/configure.in, dev/configure, dev/Makefile.am, dev/Makefile.in
#
# if you have problems, please let me know (tosa@nmr.mgh.harvard.edu)
#########################################################################
# from ------------------------------------------------
AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mri_hires_register
mri_hires_register_SOURCES=mri_hires_register.c
mri_hires_register_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_hires_register_LDFLAGS=$(OS_LDFLAGS)




