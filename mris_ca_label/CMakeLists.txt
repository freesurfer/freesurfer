project(mris_ca_label)
include_directories(${mris_ca_label_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_ca_label_SRCS
mris_ca_label.c
)


add_executable(mris_ca_label ${mris_ca_label_SRCS})
target_link_libraries(mris_ca_label ${FS_LIBS})
install(TARGETS mris_ca_label DESTINATION bin)	

