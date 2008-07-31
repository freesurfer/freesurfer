project(mri_paint)
include_directories(${mri_paint_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_paint_SRCS
mri_paint.c
)


add_executable(mri_paint ${mri_paint_SRCS})
target_link_libraries(mri_paint ${FS_LIBS})
install(TARGETS mri_paint DESTINATION bin)	

