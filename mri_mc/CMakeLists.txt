project(mri_mc)
include_directories(${mri_mc_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_mc_SRCS
build_MC_table.c
build_test_volume.c
mri_mc.c
mri_pretess.c
mris_extract_main_component.c
)


add_executable(mri_mc ${mri_mc_SRCS})
target_link_libraries(mri_mc ${FS_LIBS})
install(TARGETS mri_mc DESTINATION bin)	


