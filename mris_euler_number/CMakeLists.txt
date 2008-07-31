project(mris_euler_number)
include_directories(${mris_euler_number_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_euler_number_SRCS
mris_euler_number.c
)


add_executable(mris_euler_number ${mris_euler_number_SRCS})
target_link_libraries(mris_euler_number ${FS_LIBS})
install(TARGETS mris_euler_number DESTINATION bin)	

