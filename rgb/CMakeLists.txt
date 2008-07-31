project(rgb)
include_directories(${rgb_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(rgb_SRCS
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


add_library(rgb STATIC ${rgb_SRCS})
