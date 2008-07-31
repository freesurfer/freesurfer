project(mri_compute_seg_overlap)
include_directories(${mri_compute_seg_overlap_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_compute_seg_overlap_SRCS
mri_compute_seg_overlap.c
)


add_executable(mri_compute_seg_overlap ${mri_compute_seg_overlap_SRCS})
target_link_libraries(mri_compute_seg_overlap ${FS_LIBS})
install(TARGETS mri_compute_seg_overlap DESTINATION bin)	

