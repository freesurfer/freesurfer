include(FindPackageHandleStandardArgs)

set(MARTINOS_ANN /usr/pubsw/packages/ann/ann_1.1.1)

# find the include dir
find_path(ANN_INCLUDE_DIR PATHS ${ANN_DIR} ${MARTINOS_ANN} NAMES ANN PATH_SUFFIXES include)

# find the libraries
find_library(ANN_LIBRARIES PATHS ${ANN_DIR} ${MARTINOS_ANN} NAMES ANN PATH_SUFFIXES lib)

find_package_handle_standard_args(ANN DEFAULT_MSG ANN_LIBRARIES ANN_INCLUDE_DIR)