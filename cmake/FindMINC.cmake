# MINC Find Module

if(NOT MINC_DIR)
  set(MINC_DIR ${FS_PACKAGES_DIR}/minc/1.5)
endif()

find_path(MINC_INCLUDE_DIR HINTS ${MINC_DIR} NAMES minc.h PATH_SUFFIXES include)
find_library(VOLUME_IO_LIB HINTS ${MINC_DIR} NAMES libvolume_io.a PATH_SUFFIXES lib)
find_library(MINC_LIB HINTS ${MINC_DIR} NAMES libminc.a PATH_SUFFIXES lib)
find_package_handle_standard_args(MINC DEFAULT_MSG MINC_INCLUDE_DIR VOLUME_IO_LIB MINC_LIB)
set(MINC_LIBRARIES ${VOLUME_IO_LIB} ${MINC_LIB})
