project(jpeg)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES
  jcapimin.c
  jcapistd.c
  jctrans.c
  jcparam.c
  jdatadst.c
  jcinit.c
  jcmaster.c
  jcmarker.c
  jcmainct.c
  jcprepct.c
  jccoefct.c
  jccolor.c
  jcsample.c
  jchuff.c
  jcphuff.c
  jcdctmgr.c
  jfdctfst.c
  jfdctflt.c
  jfdctint.c
  jdapimin.c
  jdapistd.c
  jdtrans.c
  jdatasrc.c
  jdmaster.c
  jdinput.c
  jdmarker.c
  jdhuff.c
  jdphuff.c
  jdmainct.c
  jdcoefct.c
  jdpostct.c
  jddctmgr.c
  jidctfst.c
  jidctflt.c
  jidctint.c
  jidctred.c
  jdsample.c
  jdcolor.c
  jquant1.c
  jquant2.c
  jdmerge.c
  jcomapi.c
  jutils.c
  jerror.c
  jmemmgr.c
  jmemnobs.c
)

add_library(jpeg STATIC ${SOURCES})

