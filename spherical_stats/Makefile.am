AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = spherical_stats
spherical_stats_SOURCES=spherical_stats.c
spherical_stats_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
spherical_stats_LDFLAGS=$(OS_LDFLAGS)


