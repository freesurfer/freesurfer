project(expat)

add_definitions(-DHAVE_EXPAT_CONFIG_H -DBYTEORDER=${BYTEORDER})

add_library(expat STATIC xmlparse.c xmltok.c xmlrole.c)
