project(mri_normalize_tp2)
include_directories(${mri_normalize_tp2_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_normalize_tp2_SRCS
mri_normalize_tp2.c
)


add_executable(mri_normalize_tp2 ${mri_normalize_tp2_SRCS})
target_link_libraries(mri_normalize_tp2 ${FS_LIBS})
install(TARGETS mri_normalize_tp2 DESTINATION bin)	


