project(mri_partial_ribbon)
include_directories(${mri_partial_ribbon_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_partial_ribbon_SRCS
mri_partial_ribbon.c
)


add_executable(mri_partial_ribbon ${mri_partial_ribbon_SRCS})
target_link_libraries(mri_partial_ribbon ${FS_LIBS})
install(TARGETS mri_partial_ribbon DESTINATION bin)	


