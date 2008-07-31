project(mris_spharm)
include_directories(${mris_spharm_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_spharm_SRCS
mris_spharm.c
)


add_executable(mris_spharm ${mris_spharm_SRCS})
target_link_libraries(mris_spharm ${FS_LIBS})
install(TARGETS mris_spharm DESTINATION bin)	

