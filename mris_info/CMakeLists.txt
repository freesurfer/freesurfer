project(mris_info)
include_directories(${mris_info_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_info_SRCS
mris_info.cpp
)


add_executable(mris_info ${mris_info_SRCS})
target_link_libraries(mris_info ${FS_LIBS})
install(TARGETS mris_info DESTINATION bin)	

