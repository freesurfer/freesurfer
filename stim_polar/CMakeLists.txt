project(stim_polar)
include_directories(${stim_polar_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(stim_polar_SRCS
stim_polar.c
)


add_executable(stim_polar ${stim_polar_SRCS})
target_link_libraries(stim_polar ${FS_LIBS})
install(TARGETS stim_polar DESTINATION bin)	

