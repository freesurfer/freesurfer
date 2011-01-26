##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include \
        -I$(top_srcdir)/tetgen \
        $(VXL_CFLAGS)   \
        $(ITK_CFLAGS) -DHAVE_ITK_LIBS \
        $(PETSC_CFLAGS) \
        $(BOOST_CFLAGS)

AM_LDFLAGS=$(OS_LDFLAGS)

bin_PROGRAMS = dmri_spline dmri_train dmri_path_meas dmri_track dmri_trk2trk

dmri_spline_SOURCES=dmri_spline.cxx spline.cxx spline.h
dmri_spline_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
dmri_spline_LDFLAGS=$(OS_LDFLAGS)

dmri_train_SOURCES=dmri_train.cxx blood.cxx blood.h spline.cxx spline.h \
                   TrackIO.cxx TrackIO.h ByteSwap.h ErrorCode.h
dmri_train_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
dmri_train_LDFLAGS=$(OS_LDFLAGS)

dmri_path_meas_SOURCES=dmri_path_meas.cxx spline.cxx spline.h
dmri_path_meas_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
dmri_path_meas_LDFLAGS=$(OS_LDFLAGS)

dmri_track_SOURCES=dmri_track.cxx coffin.cxx coffin.h \
                   bite.cxx bite.h spline.cxx spline.h vial.cxx vial.h
dmri_track_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH)) \
                  $(top_builddir)/fem_elastic/libfem_elastic.a \
                  $(top_builddir)/tetgen/libtetgen.a 
dmri_track_LDFLAGS=$(OS_LDFLAGS)

dmri_trk2trk_SOURCES=dmri_trk2trk.cxx vial.cxx vial.h TrackIO.cxx TrackIO.h
dmri_trk2trk_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH)) \
                  $(top_builddir)/fem_elastic/libfem_elastic.a \
                  $(top_builddir)/tetgen/libtetgen.a 
dmri_trk2trk_LDFLAGS=$(OS_LDFLAGS)

EXTRA_DIST=TrackExporter.h

EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
