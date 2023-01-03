# ANN Find Module

if(NOT ANN_DIR)
  if(HOST_OS MATCHES "Ubuntu20" OR HOST_OS MATCHES "Ubuntu22")
    set(ANN_DIR ${FS_PACKAGES_DIR}/ann/1.1.2_ubuntu20)
  else()
    set(ANN_DIR ${FS_PACKAGES_DIR}/ann/1.1.2)
  endif()
endif()

find_path(ANN_INCLUDE_DIR HINTS ${ANN_DIR} NAMES ANN PATH_SUFFIXES include)
find_library(ANN_LIBRARIES HINTS ${ANN_DIR} NAMES libANN.a PATH_SUFFIXES lib)
find_package_handle_standard_args(ANN DEFAULT_MSG ANN_INCLUDE_DIR ANN_LIBRARIES)
