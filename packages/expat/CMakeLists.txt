project(expat)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_definitions(-DHAVE_EXPAT_CONFIG_H -DBYTEORDER=${BYTEORDER})

set(SOURCES
  xmlparse.c
  xmltok.c
  xmlrole.c
)

add_library(expat STATIC ${SOURCES})
