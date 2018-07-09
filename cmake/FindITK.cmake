
# if the ITK_DIR path is not supplied and we're at Martinos,
# use this default path on pubsw
set(MARTINOS_ITK_DIR /usr/pubsw/packages/itk/5.0.0)
if(NOT ITK_DIR AND EXISTS ${MARTINOS_ITK_DIR})
  set(ITK_DIR ${MARTINOS_ITK_DIR})
endif()

# if ITK_DIR was not supplied and we're building otuside of Martinos,
# download and compile ITK as an external project
if(NOT ITK_DIR)
  include(ExternalProject)
  # todo external package...
endif()

# now, actually run the find_package to set ITK variables
find_package(ITK PATHS ${ITK_DIR} NO_DEFAULT_PATH REQUIRED)
