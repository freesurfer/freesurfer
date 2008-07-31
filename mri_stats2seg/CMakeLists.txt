project(mri_stats2seg)
include_directories(${mri_stats2seg_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_stats2seg_SRCS
mri_stats2seg.c
)

SET (MRI_STATS2SEG_LIBS
${FS_LIBS}
fsgdf
)

add_executable(mri_stats2seg ${mri_stats2seg_SRCS})
target_link_libraries(mri_stats2seg ${MRI_STATS2SEG_LIBS})
install(TARGETS mri_stats2seg DESTINATION bin)	


