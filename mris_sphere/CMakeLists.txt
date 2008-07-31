project(mris_sphere)
include_directories(${mris_sphere_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_sphere_SRCS
mris_remove_negative_vertices.c
mris_sphere.c
)


add_executable(mris_sphere ${mris_sphere_SRCS})
target_link_libraries(mris_sphere ${FS_LIBS})
install(TARGETS mris_sphere DESTINATION bin)	

