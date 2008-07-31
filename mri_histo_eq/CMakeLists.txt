project(mri_histo_eq)
include_directories(${mri_histo_eq_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_histo_eq_SRCS
mri_histo_eq.c
)


add_executable(mri_histo_eq ${mri_histo_eq_SRCS})
target_link_libraries(mri_histo_eq ${FS_LIBS})
install(TARGETS mri_histo_eq DESTINATION bin)	


