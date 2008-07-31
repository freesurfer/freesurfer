project(mri_parse_sdcmdir)
include_directories(${mri_parse_sdcmdir_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_parse_sdcmdir_SRCS
mri_parse_sdcmdir.c
)


add_executable(mri_parse_sdcmdir ${mri_parse_sdcmdir_SRCS})
target_link_libraries(mri_parse_sdcmdir ${FS_LIBS})
install(TARGETS mri_parse_sdcmdir DESTINATION bin)	
