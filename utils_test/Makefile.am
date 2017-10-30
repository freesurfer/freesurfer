SUBDIRS =

AM_CFLAGS=-I$(top_srcdir)/include $(BOOST_CFLAGS)

AM_CXXFLAGS=-I$(top_srcdir)/include $(BOOST_CFLAGS) -Wall -Wextra -std=c++11

AM_LDFLAGS=

if HAVE_BOOST_LIBS

bin_PROGRAMS = Utils_Test

Utils_Test_SOURCES = boosttests.cpp\
                     mritests.cpp
Utils_Test_CXXFLAGS = -Wno-unused-parameter $(BOOST_CFLAGS) $(AM_CXXFLAGS)
Utils_Test_LDADD=  $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(BOOST_PO_LIB)
Utils_Test_LDFLAGS=$(OS_LDFLAGS)

endif
#endif HAVE_BOOST_LIBS


# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra

