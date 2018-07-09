include(ExternalProject)

# build the expat library within the freesurfer build tree
set(EXPAT_DIR ${CMAKE_BINARY_DIR}/packages/expat)
ExternalProject_Add(
  expat
  PREFIX ${EXPAT_DIR}
  URL "file://${CMAKE_SOURCE_DIR}/packages/expat-2.0.1-patch.tar.gz"
  BINARY_DIR ${EXPAT_DIR}/src/expat
  CONFIGURE_COMMAND ./buildconf.sh && ./configure --prefix=${EXPAT_DIR}
  BUILD_COMMAND make
  INSTALL_COMMAND make install
  DOWNLOAD_NO_PROGRESS true
)

# set EXPAT paths
set(EXPAT_INCLUDE_DIR ${EXPAT_DIR}/include)
set(EXPAT_LIBRARY ${EXPAT_DIR}/lib/libexpat.a)
