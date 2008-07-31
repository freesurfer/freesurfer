project(label_area)
include_directories(${label_area_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(label_area_SRCS
label_area.c
)


add_executable(label_area ${label_area_SRCS})
target_link_libraries(label_area ${FS_LIBS})
install(TARGETS label_area DESTINATION bin)	

