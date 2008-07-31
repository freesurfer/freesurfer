project(mri_fslmat_to_lta)
include_directories(${mri_fslmat_to_lta_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_fslmat_to_lta_SRCS
mri_fslmat_to_lta.c
)


add_executable(mri_fslmat_to_lta ${mri_fslmat_to_lta_SRCS})
target_link_libraries(mri_fslmat_to_lta ${FS_LIBS})
install(TARGETS mri_fslmat_to_lta DESTINATION bin)	

