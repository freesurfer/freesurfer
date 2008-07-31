project(mri_glmfit)
include_directories(${mri_glmfit_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_glmfit_SRCS
mri_glmfit.c
)

SET (mri_glmfit_LIBS
fsgdf
${FS_LIBS}
)

add_executable(mri_glmfit ${mri_glmfit_SRCS})
target_link_libraries(mri_glmfit ${mri_glmfit_LIBS})
install(TARGETS mri_glmfit DESTINATION bin)	

