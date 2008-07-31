project(mri_probedicom)
include_directories(${mri_probedicom_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_probedicom_SRCS
mri_probedicom.c
)


add_executable(mri_probedicom ${mri_probedicom_SRCS})
target_link_libraries(mri_probedicom ${FS_LIBS})
install(TARGETS mri_probedicom DESTINATION bin)	
