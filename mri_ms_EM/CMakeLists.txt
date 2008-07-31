project(mri_ms_EM)
include_directories(${mri_ms_EM_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_ms_EM_SRCS
mri_ms_EM.c
myutils.c
PoissonSolver.c
)


add_executable(mri_ms_EM ${mri_ms_EM_SRCS})
target_link_libraries(mri_ms_EM ${FS_LIBS})
install(TARGETS mri_ms_EM DESTINATION bin)	

