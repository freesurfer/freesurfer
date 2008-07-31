project(mri_interpolate)
include_directories(${mri_interpolate_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_interpolate_SRCS
mri_interpolate.c
)


add_executable(mri_interpolate ${mri_interpolate_SRCS})
target_link_libraries(mri_interpolate ${FS_LIBS})
install(TARGETS mri_interpolate DESTINATION bin)	


