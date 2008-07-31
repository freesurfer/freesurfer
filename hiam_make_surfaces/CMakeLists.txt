project(hiam_make_surfaces)
include_directories(${hiam_make_surfaces_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(hiam_make_surfaces_SRCS
hiam_make_surfaces.c
)


add_executable(hiam_make_surfaces ${hiam_make_surfaces_SRCS})
target_link_libraries(hiam_make_surfaces ${FS_LIBS})
install(TARGETS hiam_make_surfaces DESTINATION bin)	

