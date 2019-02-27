project(dicom)

if(IS_BIG_ENDIAN)
  add_definitions(-DBIG_ENDIAN_ARCHITECTURE)
else()
  add_definitions(-DLITTLE_ENDIAN_ARCHITECTURE)
endif()

add_library(dicom STATIC
  condition.c
  ctnthread.c
  dcm.c
  dcmcond.c
  dcmdict.c
  dcmsupport.c
  lst.c
)
