# Tetgen Find Module

if(NOT Tetgen_DIR)
  set(Tetgen_DIR ${FS_PACKAGES_DIR}/tetgen/1.4.1)
endif()

find_path(Tetgen_INCLUDE_DIR HINTS ${Tetgen_DIR} NAMES tetgen.h PATH_SUFFIXES include)
find_library(Tetgen_LIBRARIES HINTS ${Tetgen_DIR} NAMES libtet.a PATH_SUFFIXES lib)
find_package_handle_standard_args(Tetgen DEFAULT_MSG Tetgen_INCLUDE_DIR Tetgen_LIBRARIES)

if(EXISTS ${Tetgen_DIR}/bin/tetgen)
  install(PROGRAMS ${Tetgen_DIR}/bin/tetgen DESTINATION bin)
endif()
