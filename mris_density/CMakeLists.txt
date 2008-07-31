project(mris_density)
include_directories(${mris_density_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_density_SRCS
mris_density.c
)


add_executable(mris_density ${mris_density_SRCS})
target_link_libraries(mris_density ${FS_LIBS})
install(TARGETS mris_density DESTINATION bin)	

