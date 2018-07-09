include(ExternalProject)

# if MINC_DIR not set and at martinos:
#  set MINC_DIR to martinos

# todo set cached

if(NOT MINC_DIR)
  # if MINC_DIR is not provided, build the package within the freesurfer build tree
  # CXX=null prevents building the cpp libraries (which we don't need)
  set(MINC_DIR ${CMAKE_BINARY_DIR}/packages/minc)
  ExternalProject_Add(
    minc
    PREFIX ${MINC_DIR}
    URL "file://${CMAKE_SOURCE_DIR}/packages/minc-1.5.tar.gz"
    BINARY_DIR ${MINC_DIR}/src/minc
    CONFIGURE_COMMAND ./configure --prefix=${MINC_DIR} CPPFLAGS=-I${NETCDF_INCLUDE_DIR} LDFLAGS=-L${NETCDF_LIBRARY}
    BUILD_COMMAND make
    INSTALL_COMMAND make install
    DOWNLOAD_NO_PROGRESS true
  )
  # make sure that netcdf is built (if not already) before minc is
  add_external_dependencies(minc netcdf)
  # set minc variables
  set(MINC_FOUND)
  set(MINC_INCLUDE_DIR ${MINC_DIR}/include)
  set(MINC_LIBRARIES ${MINC_DIR}/lib/libvolume_io.a ${MINC_DIR}/lib/libminc.a)
else()
  message(FATAL_ERROR "non-local minc not yet implemented, must build within freesurfer")
endif()
