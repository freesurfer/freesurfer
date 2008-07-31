project(mri_vol2roi)
include_directories(${mri_vol2roi_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_vol2roi_SRCS
mri_vol2roi.c
)


add_executable(mri_vol2roi ${mri_vol2roi_SRCS})
target_link_libraries(mri_vol2roi ${FS_LIBS})
install(TARGETS mri_vol2roi DESTINATION bin)	


