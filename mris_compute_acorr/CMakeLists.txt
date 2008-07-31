project(mris_compute_acorr)
include_directories(${mris_compute_acorr_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_compute_acorr_SRCS
mris_compute_acorr.c
)


add_executable(mris_compute_acorr ${mris_compute_acorr_SRCS})
target_link_libraries(mris_compute_acorr ${FS_LIBS})
install(TARGETS mris_compute_acorr DESTINATION bin)	

