project(mris_aseg_distance)
include_directories(${mris_aseg_distance_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_aseg_distance_SRCS
mris_aseg_distance.c
)


add_executable(mris_aseg_distance ${mris_aseg_distance_SRCS})
target_link_libraries(mris_aseg_distance ${FS_LIBS})
install(TARGETS mris_aseg_distance DESTINATION bin)	

