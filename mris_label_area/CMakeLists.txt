project(mris_label_area)
include_directories(${mris_label_area_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_label_area_SRCS
mris_label_area.c
)


add_executable(mris_label_area ${mris_label_area_SRCS})
target_link_libraries(mris_label_area ${FS_LIBS})
install(TARGETS mris_label_area DESTINATION bin)	

