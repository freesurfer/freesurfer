project(netcdf)

include_directories(${FS_INCLUDE_DIRS})

set(SOURCES
  attr.c
  dim.c
  error.c
  libvers.c
  nc.c
  ncio.c
  ncx.c
  putget.c
  string.c
  v1hpg.c
  v2i.c
  var.c
)

add_library(netcdf STATIC ${SOURCES})
