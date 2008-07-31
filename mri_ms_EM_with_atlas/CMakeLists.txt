project(mri_ms_EM_with_atlas)
include_directories(${mri_ms_EM_with_atlas_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_ms_EM_with_atlas_SRCS
mri_ms_EM_with_atlas.c
myutils.c
PoissonSolver.c
subroutines.c
)


add_executable(mri_ms_EM_with_atlas ${mri_ms_EM_with_atlas_SRCS})
target_link_libraries(mri_ms_EM_with_atlas ${FS_LIBS})
install(TARGETS mri_ms_EM_with_atlas DESTINATION bin)	

