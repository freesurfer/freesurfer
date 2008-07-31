project(mris_thickness_diff)
include_directories(${mris_thickness_diff_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_thickness_diff_SRCS
mris_thickness_diff.c
)


add_executable(mris_thickness_diff ${mris_thickness_diff_SRCS})
target_link_libraries(mris_thickness_diff ${FS_LIBS})
install(TARGETS mris_thickness_diff DESTINATION bin)	

