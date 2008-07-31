project(mri_aparc2aseg)
include_directories(${mri_aparc2aseg_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_aparc2aseg_SRCS
mri_aparc2aseg.c
)


add_executable(mri_aparc2aseg ${mri_aparc2aseg_SRCS})
target_link_libraries(mri_aparc2aseg ${FS_LIBS})
install(TARGETS mri_aparc2aseg DESTINATION bin)	

