# Ensmallen Find Module

if(NOT ENSMALLEN_DIR)
  set(ENSMALLEN_DIR ${FS_PACKAGES_DIR}/ensmallen/1.14.0)
endif()

## FIX ME - not working
# find_path(ENSMALLEN_INCLUDE_DIRS HINTS ${ENSMALLEN_DIR} NAMES ensmallen PATH_SUFFIXES include)
set(ENSMALLEN_INCLUDE_DIRS ${ENSMALLEN_DIR}/include)
set(ENSMALLEN_INCLUDE_DIR ${ENSMALLEN_DIR}/include)

# find_library(ENSMALLEN_LIBRARIES HINTS ${ENSMALLEN_DIR} NAMES ensmallen PATH_SUFFIXES lib64 lib)
# find_package_handle_standard_args(ENSMALLEN DEFAULT_MSG ENSMALLEN_INCLUDE_DIRS ENSMALLEN_LIBRARIES)

# find_package_handle_standard_args(Ensmallen DEFAULT_MSG ENSMALLEN_INCLUDE_DIRS)
find_package_handle_standard_args(ENSMALLEN DEFAULT_MSG ENSMALLEN_INCLUDE_DIRS)

