project(mri_strip_subject_info)
include_directories(${mri_strip_subject_info_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_strip_subject_info_SRCS
mri_strip_subject_info.c
)


add_executable(mri_strip_subject_info ${mri_strip_subject_info_SRCS})
target_link_libraries(mri_strip_subject_info ${FS_LIBS})
install(TARGETS mri_strip_subject_info DESTINATION bin)	

