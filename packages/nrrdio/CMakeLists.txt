project(nrrdio)

add_definitions(-DTEEM_DIO=0 -DTEEM_QNANHIBIT=1 -DTEEM_ENDIAN=${BYTEORDER})

if(CMAKE_SIZEOF_VOID_P EQUAL 8)
  add_definitions(-DTEEM_32BIT=0)
elseif(CMAKE_SIZEOF_VOID_P EQUAL 4)
  add_definitions(-DTEEM_32BIT=1)
endif()

add_library(nrrdio STATIC
  754.c
  formatPNG.c
  accessors.c
  formatPNM.c
  array.c
  formatText.c
  arraysNrrd.c
  formatVTK.c
  axis.c
  gzio.c
  biffbiff.c
  keyvalue.c
  comment.c
  methodsNrrd.c
  defaultsNrrd.c
  miscAir.c
  dio.c
  mop.c
  encodingAscii.c
  parseAir.c
  encodingBzip2.c
  parseNrrd.c
  encoding.c
  encodingGzip.c
  encodingHex.c
  read.c
  encodingRaw.c
  reorder.c
  endianAir.c
  endianNrrd.c
  sane.c
  enum.c
  simple.c
  enumsNrrd.c
  string.c
  format.c
  subset.c
  formatEPS.c
  write.c
  formatNRRD.c
)
