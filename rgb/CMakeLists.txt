project(rgb)

include_directories(${FS_INCLUDE_DIRS})

set(SOURCES
  close.c
  filbuf.c
  flsbuf.c
  name.c
  open.c
  pix.c
  rdwr.c
  rgbutils.c
  rle.c
  row.c
)

add_library(rgb STATIC ${SOURCES})
