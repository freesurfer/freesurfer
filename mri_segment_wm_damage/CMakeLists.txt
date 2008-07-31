project(mri_segment_wm_damage)
include_directories(${mri_segment_wm_damage_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_segment_wm_damage_SRCS
mri_segment_wm_damage.c
)


add_executable(mri_segment_wm_damage ${mri_segment_wm_damage_SRCS})
target_link_libraries(mri_segment_wm_damage ${FS_LIBS})
install(TARGETS mri_segment_wm_damage DESTINATION bin)	

