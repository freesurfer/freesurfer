project(mri_ms_LDA)
include_directories(${mri_ms_LDA_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_ms_LDA_SRCS
mri_ms_LDA.c
)


add_executable(mri_ms_LDA ${mri_ms_LDA_SRCS})
target_link_libraries(mri_ms_LDA ${FS_LIBS})
install(TARGETS mri_ms_LDA DESTINATION bin)	

