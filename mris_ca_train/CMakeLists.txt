project(mris_ca_train)
include_directories(${mris_ca_train_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_ca_train_SRCS
mris_ca_train.c
)


add_executable(mris_ca_train ${mris_ca_train_SRCS})
target_link_libraries(mris_ca_train ${FS_LIBS})
install(TARGETS mris_ca_train DESTINATION bin)	
