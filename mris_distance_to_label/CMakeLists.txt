project(mris_distance_to_label)
include_directories(${mris_distance_to_label_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_distance_to_label_SRCS
mris_distance_to_label.cpp
)


add_executable(mris_distance_to_label ${mris_distance_to_label_SRCS})
target_link_libraries(mris_distance_to_label ${FS_LIBS})
install(TARGETS mris_distance_to_label DESTINATION bin)	
