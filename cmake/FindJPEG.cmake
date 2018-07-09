include(ExternalProject)

# build the jpeg library within the freesurfer build tree
set(JPEG_DIR ${CMAKE_BINARY_DIR}/packages/jpeg)
ExternalProject_Add(
  jpeg
  PREFIX ${JPEG_DIR}
  URL "file://${CMAKE_SOURCE_DIR}/packages/jpeg-6b.tar.gz"
  BINARY_DIR ${JPEG_DIR}/src/jpeg
  CONFIGURE_COMMAND ./configure --prefix=${JPEG_DIR}
  BUILD_COMMAND make
  INSTALL_COMMAND mkdir ${JPEG_DIR}/lib ${JPEG_DIR}/include && make install-lib
  DOWNLOAD_NO_PROGRESS true
)

# set JPEG paths
set(JPEG_INCLUDE_DIR ${JPEG_DIR}/include)
set(JPEG_LIBRARY ${JPEG_DIR}/lib/libjpeg.a)
