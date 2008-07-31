project(mris_convert)
include_directories(${mris_convert_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_convert_SRCS
mris_convert.c
)


add_executable(mris_convert ${mris_convert_SRCS})
target_link_libraries(mris_convert ${FS_LIBS})
install(TARGETS mris_convert DESTINATION bin)	

