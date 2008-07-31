project(mri_relabel_hypointensities)
include_directories(${mri_relabel_hypointensities_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_relabel_hypointensities_SRCS
mri_relabel_hypointensities.c
)


add_executable(mri_relabel_hypointensities ${mri_relabel_hypointensities_SRCS})
target_link_libraries(mri_relabel_hypointensities ${FS_LIBS})
install(TARGETS mri_relabel_hypointensities DESTINATION bin)	

