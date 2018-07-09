include(FindPackageHandleStandardArgs)
include(ExternalProject)

# todo if not mni_dir and @ martinos center:
set(VXL_DIR /usr/pubsw/packages/vxl/1.13.0)

if(VXL_DIR)
  # include directories
  find_path(VXL_INCLUDE_DIR PATHS ${VXL_DIR} NAMES vxl PATH_SUFFIXES include NO_DEFAULT_PATH)
  set(VXL_INCLUDE_DIRS ${VXL_INCLUDE_DIR}/vxl/core
                       ${VXL_INCLUDE_DIR}/vxl/vcl
                       ${VXL_INCLUDE_DIR}/vxl/v3p/netlib
                       ${VXL_INCLUDE_DIR}/vxl/v3p/netlib/opt)
  # find libraries
  find_library(VNL_LIBRARY PATHS ${VXL_DIR} NAMES vnl PATH_SUFFIXES lib NO_DEFAULT_PATH)
  find_library(NETLIB_LIBRARY PATHS ${VXL_DIR} NAMES netlib PATH_SUFFIXES lib NO_DEFAULT_PATH)
  find_library(V3P_NETLIB_LIBRARY PATHS ${VXL_DIR} NAMES v3p_netlib PATH_SUFFIXES lib NO_DEFAULT_PATH)
  find_library(VNL_ALGO_LIBRARY PATHS ${VXL_DIR} NAMES vnl_algo PATH_SUFFIXES lib NO_DEFAULT_PATH)
  find_library(VBL_LIBRARY PATHS ${VXL_DIR} NAMES vbl PATH_SUFFIXES lib NO_DEFAULT_PATH)
  find_library(VCL_LIBRARY PATHS ${VXL_DIR} NAMES vcl PATH_SUFFIXES lib NO_DEFAULT_PATH)
  find_library(VPL_LIBRARY PATHS ${VXL_DIR} NAMES vpl PATH_SUFFIXES lib NO_DEFAULT_PATH)
  find_library(VUL_LIBRARY PATHS ${VXL_DIR} NAMES vul PATH_SUFFIXES lib NO_DEFAULT_PATH)
  find_library(TEST_LIBRARY PATHS ${VXL_DIR} NAMES testlib PATH_SUFFIXES lib NO_DEFAULT_PATH)
  set(VXL_LIBRARIES ${VNL_LIBRARY} ${TEST_LIBRARY} ${NETLIB_LIBRARY} ${V3P_NETLIB_LIBRARY} ${VNL_ALGO_LIBRARY} ${VBL_LIBRARY} ${VCL_LIBRARY} ${VPL_LIBRARY} ${VUL_LIBRARY})
  # check package
  find_package_handle_standard_args(VXL DEFAULT_MSG VXL_LIBRARIES VXL_INCLUDE_DIR)
else()
  # build the local vxl library
endif()

# this should change because it prints inaccurate info
message(STATUS "VXL_INCLUDE_DIRS: ${VXL_INCLUDE_DIRS}")
message(STATUS "VXL_LIBRARIES: ${VXL_LIBRARIES}")
