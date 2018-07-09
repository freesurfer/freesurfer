include(FindPackageHandleStandardArgs)

set(MARTINOS_GTS /usr/pubsw/packages/gts/0.7.6)
find_path(GTS_INCLUDE_DIR PATHS ${GTS_DIR} ${MARTINOS_GTS} NAMES gts.h PATH_SUFFIXES include)
find_library(GTS_LIBRARY PATHS ${GTS_DIR} ${MARTINOS_GTS} NAMES libgts.a PATH_SUFFIXES lib)
find_package_handle_standard_args(GTS DEFAULT_MSG GTS_LIBRARY GTS_INCLUDE_DIR)
