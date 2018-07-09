project(dicom)

include_directories(
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}/include/dicom
)

SET(SOURCES
  condition.c
  ctnthread.c
  dcm.c
  dcmcond.c
  dcmdict.c
  dcmsupport.c
  lst.c
)

add_library(dicom STATIC ${SOURCES})
