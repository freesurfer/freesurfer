include(ExternalProject)

# build the jpeg library within the freesurfer build tree
set(TETGEN_DIR ${CMAKE_BINARY_DIR}/packages/tetgen)
ExternalProject_Add(
  tetgen
  PREFIX ${TETGEN_DIR}
  URL "file://${CMAKE_SOURCE_DIR}/packages/tetgen-1.4.1.tar.gz"
  BINARY_DIR ${TETGEN_DIR}/src/tetgen
  CONFIGURE_COMMAND ""
  BUILD_COMMAND make && make tetlib
  INSTALL_COMMAND cd ${TETGEN_DIR} &&
                  mkdir lib include bin &&
                  install src/tetgen/tetgen bin/ &&
                  install src/tetgen/libtet.a lib/ &&
                  install src/tetgen/tetgen.h include/
  DOWNLOAD_NO_PROGRESS true
)

# set JPEG paths
set(TETGEN_INCLUDE_DIR ${TETGEN_DIR}/include)
set(TETGEN_LIBRARY ${TETGEN_DIR}/lib/libtet.a)
