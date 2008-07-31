project(histo_compute_joint_density)
include_directories(${histo_compute_joint_density_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(histo_compute_joint_density_SRCS
histo_compute_joint_density.c
)


add_executable(histo_compute_joint_density ${histo_compute_joint_density_SRCS})
target_link_libraries(histo_compute_joint_density ${FS_LIBS})
install(TARGETS histo_compute_joint_density DESTINATION bin)	


