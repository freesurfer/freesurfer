project(mri_nlfilter)
include_directories(${mri_nlfilter_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_nlfilter_SRCS
mri_nlfilter.c
)


add_executable(mri_nlfilter ${mri_nlfilter_SRCS})
target_link_libraries(mri_nlfilter ${FS_LIBS})
install(TARGETS mri_nlfilter DESTINATION bin)	

