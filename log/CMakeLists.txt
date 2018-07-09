project(log)

include_directories(${CMAKE_SOURCE_DIR}/include)

SET(SOURCES
  congraph.c
  connect.c
  h_logz.c
  lmedian.c
  map.c
  runfuncs.c
)

add_library(log STATIC ${SOURCES})
