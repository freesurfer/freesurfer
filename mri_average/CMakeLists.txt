project(mri_average)
include_directories(${mri_average_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_average_SRCS
mri_average.c
)


add_executable(mri_average ${mri_average_SRCS})
target_link_libraries(mri_average ${FS_LIBS})
install(TARGETS mri_average DESTINATION bin)	

