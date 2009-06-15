##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mri_robust_register mri_robust_template lta_diff 
mri_robust_register_SOURCES=mri_robust_register.cpp Regression.cpp Registration.cpp CostFunctions.cpp RobustGaussian.cpp RegPowell.cpp
mri_robust_register_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_robust_register_LDFLAGS=$(OS_LDFLAGS)
mri_robust_template_SOURCES=mri_robust_template.cpp Regression.cpp Registration.cpp CostFunctions.cpp RobustGaussian.cpp
mri_robust_template_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_robust_template_LDFLAGS=$(OS_LDFLAGS)
lta_diff_SOURCES=lta_diff.cpp Regression.cpp Registration.cpp CostFunctions.cpp  RobustGaussian.cpp
lta_diff_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
lta_diff_LDFLAGS=$(OS_LDFLAGS)

EXTRA_DIST=Regression.h Registration.h Quaternion.h CostFunctions.h RobustGaussian.h RegPowell.h

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
