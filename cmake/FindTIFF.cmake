# TIFF Find Module

if(NOT TIFF_DIR)
  set(TIFF_DIR ${FS_PACKAGES_DIR}/tiff/3.6.1)
endif()

find_path(TIFF_INCLUDE_DIR HINTS ${TIFF_DIR} NAMES tiff.h PATH_SUFFIXES include)
find_library(TIFF_LIBRARIES HINTS ${TIFF_DIR} NAMES libtiff.a PATH_SUFFIXES lib)
find_package_handle_standard_args(TIFF DEFAULT_MSG TIFF_INCLUDE_DIR TIFF_LIBRARIES)
