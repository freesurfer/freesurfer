project(mris_talairach)
include_directories(${mris_talairach_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_talairach_SRCS
mris_talairach.c
)


add_executable(mris_talairach ${mris_talairach_SRCS})
target_link_libraries(mris_talairach ${FS_LIBS})
install(TARGETS mris_talairach DESTINATION bin)	


