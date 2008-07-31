project(mri_ca_normalize)
include_directories(${mri_ca_normalize_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_ca_normalize_SRCS
mri_ca_normalize.c
)


add_executable(mri_ca_normalize ${mri_ca_normalize_SRCS})
target_link_libraries(mri_ca_normalize ${FS_LIBS})
install(TARGETS mri_ca_normalize DESTINATION bin)	

