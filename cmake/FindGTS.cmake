# GTS Find Module

# GTS is only used by mris_decimate, so it's not too necessary that we supply
# it as an extra package for those that don't have it installed. Using GTS
# also requires GLIB, so we'll include that in here as well

if(NOT GTS_DIR)
  set(GTS_DIR ${FS_PACKAGES_DIR}/gts/0.7.6)
endif()

find_path(GTS_INCLUDE_DIR HINTS ${GTS_DIR} NAMES gts.h PATH_SUFFIXES include)
find_library(GTS_LIBRARY HINTS ${GTS_DIR} NAMES libgts.a PATH_SUFFIXES lib)

# let's attach the glib-2.0 directories and library as well, since they're required by GTS
find_path(GLIB_INCLUDE_DIR NAMES glib.h PATH_SUFFIXES glib-2.0)
find_path(GLIB_CONFIG_INCLUDE_DIR NAMES glibconfig.h HINTS /usr/lib64/glib-2.0 PATH_SUFFIXES include)
find_library(GLIB_LIBRARY NAMES glib-2.0)

find_package_handle_standard_args(GTS DEFAULT_MSG
  GTS_INCLUDE_DIR
  GLIB_INCLUDE_DIR
  GLIB_CONFIG_INCLUDE_DIR
  GTS_LIBRARY
  GLIB_LIBRARY
)

set(GTS_INCLUDE_DIRS ${GTS_INCLUDE_DIR} ${GLIB_INCLUDE_DIR} ${GLIB_CONFIG_INCLUDE_DIR})
set(GTS_LIBRARIES ${GTS_LIBRARY} ${GLIB_LIBRARY})
