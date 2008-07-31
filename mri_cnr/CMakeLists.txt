project(mri_cnr)
include_directories(${mri_cnr_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_cnr_SRCS
mri_cnr.c
)


add_executable(mri_cnr ${mri_cnr_SRCS})
target_link_libraries(mri_cnr ${FS_LIBS})
install(TARGETS mri_cnr DESTINATION bin)	

