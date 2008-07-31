project(mri_gdfglm)
include_directories(${mri_gdfglm_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_gdfglm_SRCS
mri_gdfglm.c
)

SET (MRI_GDFGLM_LIBS
${FS_LIBS}
fsgdf
)

add_executable(mri_gdfglm ${mri_gdfglm_SRCS})
target_link_libraries(mri_gdfglm ${MRI_GDFGLM_LIBS})
install(TARGETS mri_gdfglm DESTINATION bin)	

