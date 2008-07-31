project(mri_label_histo)
include_directories(${mri_label_histo_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_label_histo_SRCS
mri_label_histo.c
)


add_executable(mri_label_histo ${mri_label_histo_SRCS})
target_link_libraries(mri_label_histo ${FS_LIBS})
install(TARGETS mri_label_histo DESTINATION bin)	


