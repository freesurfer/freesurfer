include(ExternalProject)

# build the xml2 library within the freesurfer build tree
set(XML2_DIR ${CMAKE_BINARY_DIR}/packages/xml2)
ExternalProject_Add(
  xml2
  PREFIX ${XML2_DIR}
  URL "file://${CMAKE_SOURCE_DIR}/packages/xml2-2.7.7.tar.gz"
  BINARY_DIR ${XML2_DIR}/src/xml2
  CONFIGURE_COMMAND ./configure --prefix=${XML2_DIR} --with-sax1 --disable-shared --with-minimum
  BUILD_COMMAND make
  INSTALL_COMMAND make install
  DOWNLOAD_NO_PROGRESS true
)

# set XML2 paths
set(XML2_INCLUDE_DIR ${XML2_DIR}/include/libxml2)
set(XML2_LIBRARY ${XML2_DIR}/lib/libxml2.a)
