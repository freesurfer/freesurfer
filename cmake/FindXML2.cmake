# xml2 Find Module

if(NOT XML2_DIR)
  set(XML2_DIR ${FS_PACKAGES_DIR}/xml2/2.7.7)
endif()

find_path(XML2_INCLUDE_DIR PATHS ${XML2_DIR} NAMES libxml PATH_SUFFIXES include/libxml2)
find_library(XML2_LIBRARIES PATHS ${XML2_DIR} NAMES libxml2.a PATH_SUFFIXES lib)
find_package_handle_standard_args(XML2 DEFAULT_MSG XML2_INCLUDE_DIR XML2_LIBRARIES)
