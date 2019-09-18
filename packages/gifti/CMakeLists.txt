project(gifti)

include_directories(
  ${CMAKE_SOURCE_DIR}/packages/nifti
  ${CMAKE_SOURCE_DIR}/packages/expat
)

add_library(gifti STATIC gifti_io.c gifti_xml.c)
