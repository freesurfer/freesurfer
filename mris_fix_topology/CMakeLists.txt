project(mris_fix_topology)
include_directories(${mris_fix_topology_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_fix_topology_SRCS
mris_fix_topology.c
mris_remove_intersection.c
)


add_executable(mris_fix_topology ${mris_fix_topology_SRCS})
target_link_libraries(mris_fix_topology ${FS_LIBS})
install(TARGETS mris_fix_topology DESTINATION bin)	

