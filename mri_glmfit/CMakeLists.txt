project(mri_glmfit)
include_directories(${mri_glmfit_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_glmfit_SRCS
mri_glmfit.c
${mri_glmfit_SOURCE_DIR}/../fsgdf/fsgdf.c
)

SET (MRI_GLMFIT_LIBS
${FS_LIBS}
fsgdf
)

add_executable(mri_glmfit ${mri_glmfit_SRCS})
target_link_libraries(mri_glmfit ${MRI_GLMFIT_LIBS})
install(TARGETS mri_glmfit DESTINATION bin)	

