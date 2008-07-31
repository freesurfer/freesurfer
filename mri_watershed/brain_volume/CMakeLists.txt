project(brain_volume)
include_directories(${brain_volume_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(brain_volume_SRCS
mri_brain_volume.cpp
)


add_executable(brain_volume ${brain_volume_SRCS})
target_link_libraries(brain_volume ${FS_LIBS})
install(TARGETS brain_volume DESTINATION bin)	

add_subdirectory(test)
