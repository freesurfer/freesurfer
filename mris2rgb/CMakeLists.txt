project(mris2rgb)
include_directories(${mris2rgb_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris2rgb_SRCS
mris2rgb.c
)


add_executable(mris2rgb ${mris2rgb_SRCS})
target_link_libraries(mris2rgb ${FS_LIBS})
install(TARGETS mris2rgb DESTINATION bin)	

