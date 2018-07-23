# Expat Find Module

if(NOT Expat_DIR)
  set(Expat_DIR ${FS_PACKAGES_DIR}/expat/2.0.1)
endif()

find_path(Expat_INCLUDE_DIR HINTS ${Expat_DIR} NAMES expat.h PATH_SUFFIXES include)
find_library(Expat_LIBRARIES HINTS ${Expat_DIR} NAMES libexpat.a PATH_SUFFIXES lib)
find_package_handle_standard_args(Expat DEFAULT_MSG Expat_INCLUDE_DIR Expat_LIBRARIES)
