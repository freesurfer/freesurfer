project(mri_strip_nonwhite)
include_directories(${mri_strip_nonwhite_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_strip_nonwhite_SRCS
mri_strip_nonwhite.c
)


add_executable(mri_strip_nonwhite ${mri_strip_nonwhite_SRCS})
target_link_libraries(mri_strip_nonwhite ${FS_LIBS})
install(TARGETS mri_strip_nonwhite DESTINATION bin)	

