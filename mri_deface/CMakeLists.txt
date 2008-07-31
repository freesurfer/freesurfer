project(mri_deface)
include_directories(${mri_deface_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_deface_SRCS
mri_deface.c
)


add_executable(mri_deface ${mri_deface_SRCS})
target_link_libraries(mri_deface ${FS_LIBS})
install(TARGETS mri_deface DESTINATION bin)	
