## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=$(GLUT_CFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = test_window_env
test_window_env_SOURCES=test_window_env.c
test_window_env_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH)) $(GLUT_LIBS)  
test_window_env_LDFLAGS=

