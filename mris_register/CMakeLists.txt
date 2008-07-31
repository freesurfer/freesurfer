project(mris_register)
include_directories(${mris_register_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_register_SRCS
mris_register.c
)


add_executable(mris_register ${mris_register_SRCS})
target_link_libraries(mris_register ${FS_LIBS})
install(TARGETS mris_register DESTINATION bin)	


