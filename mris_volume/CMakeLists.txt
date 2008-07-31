project(mris_volume)
include_directories(${mris_volume_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_volume_SRCS
mris_volume.c
mris_wm_volume.c
)


add_executable(mris_volume ${mris_volume_SRCS})
target_link_libraries(mris_volume ${FS_LIBS})
install(TARGETS mris_volume DESTINATION bin)	


