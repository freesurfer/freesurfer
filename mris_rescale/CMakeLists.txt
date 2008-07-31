project(mris_rescale)
include_directories(${mris_rescale_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_rescale_SRCS
mris_rescale.c
)


add_executable(mris_rescale ${mris_rescale_SRCS})
target_link_libraries(mris_rescale ${FS_LIBS})
install(TARGETS mris_rescale DESTINATION bin)	

