project(mris_compute_optimal_kernel)
include_directories(${mris_compute_optimal_kernel_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_compute_optimal_kernel_SRCS
mris_compute_optimal_kernel.c
)


add_executable(mris_compute_optimal_kernel ${mris_compute_optimal_kernel_SRCS})
target_link_libraries(mris_compute_optimal_kernel ${FS_LIBS})
install(TARGETS mris_compute_optimal_kernel DESTINATION bin)	

