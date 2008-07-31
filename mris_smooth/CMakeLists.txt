project(mris_smooth)
include_directories(${mris_smooth_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_smooth_SRCS
mris_smooth.c
)


add_executable(mris_smooth ${mris_smooth_SRCS})
target_link_libraries(mris_smooth ${FS_LIBS})
install(TARGETS mris_smooth DESTINATION bin)	

