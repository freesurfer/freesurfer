project(mri_morphology)
include_directories(${mri_morphology_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_morphology_SRCS
mri_morphology.c
)


add_executable(mri_morphology ${mri_morphology_SRCS})
target_link_libraries(mri_morphology ${FS_LIBS})
install(TARGETS mri_morphology DESTINATION bin)	

