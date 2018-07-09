include(FindPackageHandleStandardArgs)

find_path(GLIB2_INCLUDE NAMES glib.h PATH_SUFFIXES glib-2.0)
find_path(GLIB2_CONFIG NAMES glibconfig.h PATH_SUFFIXES glib-2.0/include)
find_library(GLIB2_LIBRARY NAMES libglib-2.0.a)
find_package_handle_standard_args(GLIB2 DEFAULT_MSG GLIB2_LIBRARY GLIB2_INCLUDE GLIB2_CONFIG)
set(GLIB2_INCLUDE_DIRS ${GLIB2_INCLUDE} ${GLIB2_CONFIG})
