##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include $(VTK_CXXFLAGS) -Wno-deprecated

if HAVE_VTK_LIBS
    bin_PROGRAMS = mris_decimate mris_decimate_gui
    mris_decimate_SOURCES=mris_decimate.cpp mris_decimate.h
    mris_decimate_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(VTK_LIBS) -lvtkHybrid -lvtkWidgets \
        -lvtkGraphics \
        -lvtkCommon  \
        -lvtkImaging \
        -lvtkFiltering
    mris_decimate_LDFLAGS=$(OS_LDFLAGS)

    mris_decimate_gui_SOURCES=mris_decimate.cpp mris_decimate.h
    mris_decimate_gui_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(VTK_LIBS) -lvtkHybrid -lvtkWidgets \
        -lvtkGraphics \
        -lvtkCommon  \
        -lvtkImaging \
        -lvtkFiltering \
        -lvtkRendering \
        -lvtkIO
    mris_decimate_gui_LDFLAGS=$(OS_LDFLAGS)
    mris_decimate_gui_CXXFLAGS=$(AM_CXXFLAGS) -DBUILT_IN_VISUALIZER

endif

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
