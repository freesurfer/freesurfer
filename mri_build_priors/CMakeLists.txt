project(mri_build_priors)
include_directories(${mri_build_priors_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_build_priors_SRCS
mri_build_priors.c
)


add_executable(mri_build_priors ${mri_build_priors_SRCS})
target_link_libraries(mri_build_priors ${FS_LIBS})
install(TARGETS mri_build_priors DESTINATION bin)	

