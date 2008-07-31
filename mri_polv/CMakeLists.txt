project(mri_polv)
include_directories(${mri_polv_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_polv_SRCS
mri_polv.c
)


add_executable(mri_polv ${mri_polv_SRCS})
target_link_libraries(mri_polv ${FS_LIBS})
install(TARGETS mri_polv DESTINATION bin)	

