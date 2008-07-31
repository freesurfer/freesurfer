project(mri_train)
include_directories(${mri_train_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_train_SRCS
mri_train.c
)


add_executable(mri_train ${mri_train_SRCS})
target_link_libraries(mri_train ${FS_LIBS})
install(TARGETS mri_train DESTINATION bin)	

