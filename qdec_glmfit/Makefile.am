AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-Wno-deprecated -I$(top_srcdir)/include \
	-I$(top_srcdir)/qdecproject

bin_PROGRAMS = qdec_glmfit
qdec_glmfit_SOURCES=qdec_glmfit.cxx
qdec_glmfit_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH)) \
	$(top_srcdir)/qdecproject/libqdecproject.a
qdec_glmfit_LDFLAGS=$(OS_LDFLAGS)

EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
