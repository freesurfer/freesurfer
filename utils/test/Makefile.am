#
# Makefile.am
#
AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

check_PROGRAMS=mnitest i2rtest icotest surftest extest mghxform inftest difftool mri_copy_params mri_remove checkanalyze surf2surf
TESTS=mnitest i2rtest icotest extest mghxform inftest checkanalyze

mnitest_SOURCES=mnitest.cpp
i2rtest_SOURCES=i2rtest.cpp
icotest_SOURCES=icotest.cpp
surftest_SOURCES=surftest.cpp
extest_SOURCES=extest.cpp
mghxform_SOURCES=mghxform.cpp
inftest_SOURCES=inftest.cpp
difftool_SOURCES=difftool.cpp
mri_copy_params_SOURCES=mri_copy_params.cpp
mri_remove_SOURCES=mri_remove.cpp
checkanalyze_SOURCES=checkanalyze.cpp
surf2surf_SOURCES=surf2surf.cpp


LDADD=$(addprefix $(top_srcdir)/, $(LIBS_MGH))

