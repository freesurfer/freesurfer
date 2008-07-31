project(mris_jacobian)
include_directories(${mris_jacobian_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_jacobian_SRCS
mris_jacobian.c
)


add_executable(mris_jacobian ${mris_jacobian_SRCS})
target_link_libraries(mris_jacobian ${FS_LIBS})
install(TARGETS mris_jacobian DESTINATION bin)	

