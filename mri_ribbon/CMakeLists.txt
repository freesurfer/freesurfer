project(mri_ribbon)
include_directories(${mri_ribbon_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_ribbon_SRCS
mri_ribbon.c
)


add_executable(mri_ribbon ${mri_ribbon_SRCS})
target_link_libraries(mri_ribbon ${FS_LIBS})
install(TARGETS mri_ribbon DESTINATION bin)	

