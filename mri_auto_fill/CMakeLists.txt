project(mri_auto_fill)
include_directories(${mri_auto_fill_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_auto_fill_SRCS
mri_auto_fill.c
)


add_executable(mri_auto_fill ${mri_auto_fill_SRCS})
target_link_libraries(mri_auto_fill ${FS_LIBS})
install(TARGETS mri_auto_fill DESTINATION bin)	

