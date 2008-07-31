project(mri_wmfilter)
include_directories(${mri_wmfilter_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_wmfilter_SRCS
mri_wmfilter.c
MRIio.c
)


add_executable(mri_wmfilter ${mri_wmfilter_SRCS})
target_link_libraries(mri_wmfilter ${FS_LIBS})
install(TARGETS mri_wmfilter DESTINATION bin)	

add_subdirectory(old)
