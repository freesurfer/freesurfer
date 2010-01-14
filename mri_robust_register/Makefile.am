##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include $(VXL_CFLAGS)

bin_PROGRAMS = mri_robust_register mri_robust_template lta_diff mri_create_tests mri_gradient_info
mri_robust_register_SOURCES=Regression.cpp \
	CostFunctions.cpp RobustGaussian.cpp MyMatrix.cpp MyMRI.cpp \
	Quaternion.cpp Registration.cpp RegPowell.cpp mri_robust_register.cpp
mri_robust_register_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_robust_register_LDFLAGS=$(OS_LDFLAGS)
mri_robust_template_SOURCES=Regression.cpp Registration.cpp \
	CostFunctions.cpp RobustGaussian.cpp MyMatrix.cpp MyMRI.cpp Quaternion.cpp \
	MultiRegistration.cpp mri_robust_template.cpp 
mri_robust_template_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_robust_template_LDFLAGS=$(OS_LDFLAGS)
lta_diff_SOURCES=lta_diff.cpp Regression.cpp Registration.cpp \
	CostFunctions.cpp  RobustGaussian.cpp MyMatrix.cpp MyMRI.cpp Quaternion.cpp
lta_diff_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
lta_diff_LDFLAGS=$(OS_LDFLAGS)
mri_create_tests_SOURCES=mri_create_tests.cpp Quaternion.cpp MyMatrix.cpp MyMRI.cpp
mri_create_tests_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_create_tests_LDFLAGS=$(OS_LDFLAGS)
mri_gradient_info_SOURCES=mri_gradient_info.cpp
mri_gradient_info_LDADD = $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_gradient_info_LDFLAGS=$(OS_LDFLAGS)

EXTRA_DIST=Regression.h Registration.h MultiRegistration.h Quaternion.h \
	CostFunctions.h RobustGaussian.h RegPowell.h MyMatrix.h MyMRI.h

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
