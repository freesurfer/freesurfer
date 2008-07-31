project(mri_aparc2wmseg)
include_directories(${mri_aparc2wmseg_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_aparc2wmseg_SRCS
mri_aparc2wmseg.c
)


add_executable(mri_aparc2wmseg ${mri_aparc2wmseg_SRCS})
target_link_libraries(mri_aparc2wmseg ${FS_LIBS})
install(TARGETS mri_aparc2wmseg DESTINATION bin)	

