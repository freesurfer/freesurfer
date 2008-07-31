project(mri_volsynth)
include_directories(${mri_volsynth_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_volsynth_SRCS
mri_volsynth.c
)


add_executable(mri_volsynth ${mri_volsynth_SRCS})
target_link_libraries(mri_volsynth ${FS_LIBS})
install(TARGETS mri_volsynth DESTINATION bin)	

