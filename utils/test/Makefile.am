#
# Makefile.am
#
AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

check_PROGRAMS=mnitest i2rtest icotest surftest extest mghxform inftest difftool tixtest mri_copy_params mri_remove mattest mri_resample rotate_label_surf LTAinvert checkanalyze
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
mri_resample_SOURCES=mri_resample.cpp
mri_resample_LDFLAGS=$(OS_LDFLAGS)
mattest_SOURCES=mattest.c
tixtest_SOURCES=tixtest.c
rotate_label_surf_SOURCES=rotate_label_surf.cpp
LTAinvert_SOURCES=LTAinvert.cpp
checkanalyze_SOURCES=checkanalyze.cpp

tixtest_LDFLAGS = $(LDFLAGS_TCL)
tixtest_LDADD=$(LIBS_TCL)

LDADD=$(addprefix $(top_srcdir)/, $(LIBS_MGH))

