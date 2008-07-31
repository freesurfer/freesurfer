project(mri_twoclass)
include_directories(${mri_twoclass_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_twoclass_SRCS
mri_twoclass.c
)


add_executable(mri_twoclass ${mri_twoclass_SRCS})
target_link_libraries(mri_twoclass ${FS_LIBS})
install(TARGETS mri_twoclass DESTINATION bin)	

