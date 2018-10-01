# Jpeg Find Module

if(NOT JPEG_DIR)
  set(JPEG_DIR ${FS_PACKAGES_DIR}/jpeg/6b)
endif()

find_path(JPEG_INCLUDE_DIR HINTS ${JPEG_DIR} NAMES jpeglib.h PATH_SUFFIXES include)
find_library(JPEG_LIBRARIES HINTS ${JPEG_DIR} NAMES libjpeg.a PATH_SUFFIXES lib)
find_package_handle_standard_args(JPEG DEFAULT_MSG JPEG_INCLUDE_DIR JPEG_LIBRARIES)
