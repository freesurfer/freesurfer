project(mri_normalize)
include_directories(${mri_normalize_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_normalize_SRCS
mri_normalize.c
)


add_executable(mri_normalize ${mri_normalize_SRCS})
target_link_libraries(mri_normalize ${FS_LIBS})
install(TARGETS mri_normalize DESTINATION bin)	

