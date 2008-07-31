project(mri_cc)
include_directories(${mri_cc_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_cc_SRCS
mri_cc.c
)


add_executable(mri_cc ${mri_cc_SRCS})
target_link_libraries(mri_cc ${FS_LIBS})
install(TARGETS mri_cc DESTINATION bin)	

