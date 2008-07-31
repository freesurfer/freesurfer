project(mri_segstats)
include_directories(${mri_segstats_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_segstats_SRCS
mri_segstats.c
)


add_executable(mri_segstats ${mri_segstats_SRCS})
target_link_libraries(mri_segstats ${FS_LIBS})
install(TARGETS mri_segstats DESTINATION bin)	

