# GTS Find Module
# note: only used by mris_decimate

if(NOT GTS_DIR)
  set(GTS_DIR ${FS_PACKAGES_DIR}/gts/0.7.6)
endif()

find_path(GTS_INCLUDE_DIR HINTS ${GTS_DIR} NAMES gts.h PATH_SUFFIXES include)
find_library(GTS_LIBRARY HINTS ${GTS_DIR} NAMES libgts.a PATH_SUFFIXES lib)

# Let's include glib-2.0 here since it's required by GTS.
# Correctly linking glib on mac is very messy, so all the glib
# headers and libraries (and the libintl.a dependency) have been copied
# into the GTS installation

find_path(GLIB_INCLUDE_DIR NAMES glib.h HINTS ${GTS_DIR} PATH_SUFFIXES glib-2.0)
find_path(GLIB_CONFIG_INCLUDE_DIR NAMES glibconfig.h HINTS ${GTS_DIR}/include /usr/lib64 PATH_SUFFIXES glib-2.0/include)
find_library(GLIB_LIBRARY NAMES libglib-2.0.a glib-2.0 HINTS ${GTS_DIR}/lib)

if(APPLE)
  find_library(INTL_LIBRARY NAMES libintl.a HINTS ${GTS_DIR} PATH_SUFFIXES lib)
  find_library(ICONV_LIBRARY NAMES iconv)
  set(APPLE_GLIB_DEPENDENCIES ICONV_LIBRARY INTL_LIBRARY)
endif()

find_package_handle_standard_args(GTS DEFAULT_MSG
  GTS_INCLUDE_DIR
  GTS_LIBRARY
  GLIB_INCLUDE_DIR
  GLIB_CONFIG_INCLUDE_DIR
  GLIB_LIBRARY
  ${APPLE_GLIB_DEPENDENCIES}
)

set(GTS_INCLUDE_DIRS ${GTS_INCLUDE_DIR} ${GLIB_INCLUDE_DIR} ${GLIB_CONFIG_INCLUDE_DIR})
set(GTS_LIBRARIES ${GTS_LIBRARY} ${GLIB_LIBRARY} ${ICONV_LIBRARY} ${INTL_LIBRARY})
