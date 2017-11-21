project(mri_exclude_minc_tests)
include_directories(${mri_exclude_minc_tests_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_exclude_minc_tests_SRCS
mri_exclude_minc_tests.c
)


add_executable(mri_exclude_minc_tests ${mri_exclude_minc_tests_SRCS})
target_link_libraries(mri_exclude_minc_tests ${FS_LIBS})
install(TARGETS mri_exclude_minc_tests DESTINATION bin)	
