project(mri_extract_conditions)
include_directories(${mri_extract_conditions_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_extract_conditions_SRCS
mri_extract_conditions.c
)


add_executable(mri_extract_conditions ${mri_extract_conditions_SRCS})
target_link_libraries(mri_extract_conditions ${FS_LIBS})
install(TARGETS mri_extract_conditions DESTINATION bin)	

