include(ExternalProject)

# build the tiff library within the freesurfer build tree
set(TIFF_DIR ${CMAKE_BINARY_DIR}/packages/tiff)
ExternalProject_Add(
  tiff
  PREFIX ${TIFF_DIR}
  URL "file://${CMAKE_SOURCE_DIR}/packages/tiff-3.6.1.tar.gz"
  BINARY_DIR ${TIFF_DIR}/src/tiff
  CONFIGURE_COMMAND ./configure --prefix=${TIFF_DIR} --noninteractive --with-PARAM=DSO=no
  BUILD_COMMAND make
  INSTALL_COMMAND make install
  DOWNLOAD_NO_PROGRESS true
)

# set TIFF paths
set(TIFF_INCLUDE_DIR ${TIFF_DIR}/include)
set(TIFF_LIBRARY ${TIFF_DIR}/lib/libtiff.a)
