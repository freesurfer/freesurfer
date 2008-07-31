project(mri_ca_train)
include_directories(${mri_ca_train_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_ca_train_SRCS
mri_ca_train.c
)


add_executable(mri_ca_train ${mri_ca_train_SRCS})
target_link_libraries(mri_ca_train ${FS_LIBS})
install(TARGETS mri_ca_train DESTINATION bin)	

