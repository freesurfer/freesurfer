project(mris_sample_label)
include_directories(${mris_sample_label_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_sample_label_SRCS
mris_sample_label.c
)


add_executable(mris_sample_label ${mris_sample_label_SRCS})
target_link_libraries(mris_sample_label ${FS_LIBS})
install(TARGETS mris_sample_label DESTINATION bin)	

