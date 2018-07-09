include(ExternalProject)

# if NETCDF_DIR not set and at martinos:
#  set NETCDF_DIR to martinos

# todo set cached

if(NOT NETCDF_DIR)
  # if NETCDF_DIR is not provided, build the package within the freesurfer build tree
  # CXX=null prevents building the cpp libraries (which we don't need)
  set(NETCDF_DIR ${CMAKE_BINARY_DIR}/packages/netcdf)
  ExternalProject_Add(
    netcdf
    PREFIX ${NETCDF_DIR}
    URL "file://${CMAKE_SOURCE_DIR}/packages/netcdf-3.6.0-p1.tar.gz"
    BINARY_DIR ${NETCDF_DIR}/src/netcdf/src
    CONFIGURE_COMMAND ./configure --prefix=${NETCDF_DIR} CXX=null
    BUILD_COMMAND make
    INSTALL_COMMAND make install
    DOWNLOAD_NO_PROGRESS true
  )
  # set netcdf variables
  set(NETCDF_FOUND)
  set(NETCDF_INCLUDE_DIR ${NETCDF_DIR}/include)
  set(NETCDF_LIBRARY ${NETCDF_DIR}/lib/libnetcdf.a)
else()
  message(FATAL_ERROR "non-local netcdf not yet implemented, must build within freesurfer")
endif()
