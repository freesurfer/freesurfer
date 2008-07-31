project(mri_mcsim)
include_directories(${mri_mcsim_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_mcsim_SRCS
mri_mcsim.c
)


add_executable(mri_mcsim ${mri_mcsim_SRCS})
target_link_libraries(mri_mcsim ${FS_LIBS})
install(TARGETS mri_mcsim DESTINATION bin)	

