project(mris_spherical_average)
include_directories(${mris_spherical_average_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_spherical_average_SRCS
mris_spherical_average.c
)


add_executable(mris_spherical_average ${mris_spherical_average_SRCS})
target_link_libraries(mris_spherical_average ${FS_LIBS})
install(TARGETS mris_spherical_average DESTINATION bin)	

