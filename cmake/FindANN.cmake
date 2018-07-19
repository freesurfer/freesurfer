# ANN Find Module

if(NOT ANN_DIR)
  set(ANN_DIR ${FS_PACKAGES_DIR}/ann/ann_1.1.1)
endif()

find_path(ANN_INCLUDE_DIR PATHS ${ANN_DIR} NAMES ANN PATH_SUFFIXES include)
find_library(ANN_LIBRARIES PATHS ${ANN_DIR} NAMES libANN.a PATH_SUFFIXES lib)
find_package_handle_standard_args(ANN DEFAULT_MSG ANN_INCLUDE_DIR ANN_LIBRARIES)
