project(mri_evaluate_morph)
include_directories(${mri_evaluate_morph_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_evaluate_morph_SRCS
mri_evaluate_morph.c
)


add_executable(mri_evaluate_morph ${mri_evaluate_morph_SRCS})
target_link_libraries(mri_evaluate_morph ${FS_LIBS})
install(TARGETS mri_evaluate_morph DESTINATION bin)	
