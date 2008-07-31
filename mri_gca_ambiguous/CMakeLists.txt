project(mri_gca_ambiguous)
include_directories(${mri_gca_ambiguous_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_gca_ambiguous_SRCS
mri_gca_ambiguous.c
)


add_executable(mri_gca_ambiguous ${mri_gca_ambiguous_SRCS})
target_link_libraries(mri_gca_ambiguous ${FS_LIBS})
install(TARGETS mri_gca_ambiguous DESTINATION bin)	


