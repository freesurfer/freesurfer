project(mri_extract_largest_CC)
include_directories(${mri_extract_largest_CC_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_extract_largest_CC_SRCS
mri_extract_largest_CC.c
myutils.c
subroutines.c
)


add_executable(mri_extract_largest_CC ${mri_extract_largest_CC_SRCS})
target_link_libraries(mri_extract_largest_CC ${FS_LIBS})
install(TARGETS mri_extract_largest_CC DESTINATION bin)	

