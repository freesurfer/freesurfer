project(mri_voldiff)
include_directories(${mri_voldiff_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_voldiff_SRCS
mri_voldiff.c
)


add_executable(mri_voldiff ${mri_voldiff_SRCS})
target_link_libraries(mri_voldiff ${FS_LIBS})
install(TARGETS mri_voldiff DESTINATION bin)	
