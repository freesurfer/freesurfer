project(netcdf_3_6_0_p1)

include_directories(
  ${FS_INCLUDE_DIRS}
)

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

add_library(netcdf_3_6_0_p1 STATIC ${SOURCES})
