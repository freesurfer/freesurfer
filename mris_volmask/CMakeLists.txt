project(mris_volmask)
include_directories(${mris_volmask_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_volmask_SRCS
cmd_line_interface.cpp
mris_volmask.c
mris_volmask.cpp
)


add_executable(mris_volmask ${mris_volmask_SRCS})
target_link_libraries(mris_volmask ${FS_LIBS})
install(TARGETS mris_volmask DESTINATION bin)	

