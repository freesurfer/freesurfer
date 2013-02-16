##
## Makefile.am 
##

if HAVE_ITK_LIBS
if HAVE_PETSC_LIBS
if HAVE_BOOST_LIBS

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include \
        -I$(top_srcdir)/fem_elastic \
        -I$(top_srcdir)/tetgen \
        $(VXL_CFLAGS)   \
        $(ITK_CFLAGS) -DHAVE_ITK_LIBS \
        $(PETSC_CFLAGS) \
        $(BOOST_CFLAGS)

AM_LDFLAGS=$(OS_LDFLAGS)

bin_PROGRAMS = dmri_spline dmri_trk2trk dmri_vox2vox \
               dmri_train dmri_paths dmri_pathstats dmri_mergepaths

dmri_spline_SOURCES=dmri_spline.cxx spline.cxx spline.h
dmri_spline_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
dmri_spline_LDFLAGS=$(OS_LDFLAGS)

dmri_trk2trk_SOURCES=dmri_trk2trk.cxx vial.cxx vial.h TrackIO.cxx TrackIO.h
dmri_trk2trk_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH)) \
                  $(top_builddir)/fem_elastic/libfem_elastic.a \
                  $(top_builddir)/tetgen/libtetgen.a 
dmri_trk2trk_LDFLAGS=$(OS_LDFLAGS)

dmri_vox2vox_SOURCES=dmri_vox2vox.cxx vial.cxx vial.h
dmri_vox2vox_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH)) \
                  $(top_builddir)/fem_elastic/libfem_elastic.a \
                  $(top_builddir)/tetgen/libtetgen.a 
dmri_vox2vox_LDFLAGS=$(OS_LDFLAGS)

dmri_train_SOURCES=dmri_train.cxx blood.cxx blood.h spline.cxx spline.h \
                   vial.cxx vial.h \
                   TrackIO.cxx TrackIO.h ByteSwap.h ErrorCode.h
dmri_train_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH)) \
                  $(top_builddir)/fem_elastic/libfem_elastic.a \
                  $(top_builddir)/tetgen/libtetgen.a 
dmri_train_LDFLAGS=$(OS_LDFLAGS)

dmri_paths_SOURCES=dmri_paths.cxx coffin.cxx coffin.h \
                   bite.cxx bite.h spline.cxx spline.h vial.cxx vial.h
dmri_paths_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH)) \
                  $(top_builddir)/fem_elastic/libfem_elastic.a \
                  $(top_builddir)/tetgen/libtetgen.a 
dmri_paths_LDFLAGS=$(OS_LDFLAGS)

dmri_pathstats_SOURCES=dmri_pathstats.cxx spline.cxx spline.h \
                       blood.cxx blood.h \
                       vial.cxx vial.h \
                       TrackIO.cxx TrackIO.h ByteSwap.h ErrorCode.h
dmri_pathstats_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH)) \
                      $(top_builddir)/fem_elastic/libfem_elastic.a \
                      $(top_builddir)/tetgen/libtetgen.a 
dmri_pathstats_LDFLAGS=$(OS_LDFLAGS)

dmri_mergepaths_SOURCES=dmri_mergepaths.cxx
dmri_mergepaths_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
dmri_mergepaths_LDFLAGS=$(OS_LDFLAGS)

# endif HAVE_ITK_LIBS, HAVE_PETSC_LIBS, HAVE_BOOST_LIBS
endif
endif
endif

EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
