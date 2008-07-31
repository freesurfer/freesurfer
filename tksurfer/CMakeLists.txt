project(tksurfer)
include_directories(${tksurfer_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(tksurfer_SRCS
tksurfer.c
)


add_executable(tksurfer ${tksurfer_SRCS})
target_link_libraries(tksurfer ${FS_LIBS})
install(TARGETS tksurfer DESTINATION bin)	

add_subdirectory(test)
