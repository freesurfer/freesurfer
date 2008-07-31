project(nmovie)
include_directories(${nmovie_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(nmovie_SRCS
nmovie.c
)


add_executable(nmovie ${nmovie_SRCS})
target_link_libraries(nmovie ${FS_LIBS})
install(TARGETS nmovie DESTINATION bin)	

