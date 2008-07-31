project(mri_fill)
include_directories(${mri_fill_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_fill_SRCS
mri_fill.c
myutil.c
subroutines.c
)


add_executable(mri_fill ${mri_fill_SRCS})
target_link_libraries(mri_fill ${FS_LIBS})
install(TARGETS mri_fill DESTINATION bin)	

