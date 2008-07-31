project(unix)
include_directories(${unix_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(unix_SRCS
cuserid.c
ftime.c
glStub.c
oglutil.c
texfont.c
unix_thread.c
windiag.c
xvmri.c
xvutil.c
xwin.c
xwindow.c
)


add_library(unix STATIC ${unix_SRCS})
