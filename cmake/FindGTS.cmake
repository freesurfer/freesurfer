# GTS Find Module
# note: only used by mris_decimate

if(NOT GTS_DIR)
  set(GTS_DIR ${FS_PACKAGES_DIR}/gts/0.7.6)
endif()

find_path(GTS_INCLUDE_DIR HINTS ${GTS_DIR} NAMES gts.h PATH_SUFFIXES include)
find_library(GTS_LIBRARY HINTS ${GTS_DIR} NAMES libgts.a PATH_SUFFIXES lib)

# let's include glib-2.0 here as well since it's required by GTS
find_path(GLIB_INCLUDE_DIR NAMES glib.h PATH_SUFFIXES glib-2.0)
find_path(GLIB_CONFIG_INCLUDE_DIR NAMES glibconfig.h HINTS /usr/lib64 /usr/local/lib PATH_SUFFIXES glib-2.0/include)
find_library(GLIB_LIBRARY NAMES libglib-2.0.a glib-2.0)

find_package_handle_standard_args(GTS DEFAULT_MSG
  GTS_INCLUDE_DIR
  GTS_LIBRARY
  GLIB_INCLUDE_DIR
  GLIB_CONFIG_INCLUDE_DIR
  GLIB_LIBRARY
)

set(GTS_INCLUDE_DIRS ${GTS_INCLUDE_DIR} ${GLIB_INCLUDE_DIR} ${GLIB_CONFIG_INCLUDE_DIR})
set(GTS_LIBRARIES ${GTS_LIBRARY} ${GLIB_LIBRARY})
