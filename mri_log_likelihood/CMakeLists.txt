project(mri_log_likelihood)
include_directories(${mri_log_likelihood_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_log_likelihood_SRCS
mri_log_likelihood.c
)


add_executable(mri_log_likelihood ${mri_log_likelihood_SRCS})
target_link_libraries(mri_log_likelihood ${FS_LIBS})
install(TARGETS mri_log_likelihood DESTINATION bin)	

