project(mri_remove_neck)
include_directories(${mri_remove_neck_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_remove_neck_SRCS
mri_remove_neck.c
)


add_executable(mri_remove_neck ${mri_remove_neck_SRCS})
target_link_libraries(mri_remove_neck ${FS_LIBS})
install(TARGETS mri_remove_neck DESTINATION bin)	


