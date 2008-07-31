project(mris_label_mode)
include_directories(${mris_label_mode_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_label_mode_SRCS
mris_label_mode.c
)


add_executable(mris_label_mode ${mris_label_mode_SRCS})
target_link_libraries(mris_label_mode ${FS_LIBS})
install(TARGETS mris_label_mode DESTINATION bin)	

