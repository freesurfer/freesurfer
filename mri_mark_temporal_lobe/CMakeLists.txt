project(mri_mark_temporal_lobe)
include_directories(${mri_mark_temporal_lobe_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_mark_temporal_lobe_SRCS
mri_mark_temporal_lobe.c
)


add_executable(mri_mark_temporal_lobe ${mri_mark_temporal_lobe_SRCS})
target_link_libraries(mri_mark_temporal_lobe ${FS_LIBS})
install(TARGETS mri_mark_temporal_lobe DESTINATION bin)	

