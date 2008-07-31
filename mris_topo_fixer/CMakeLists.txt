project(mris_topo_fixer)
include_directories(${mris_topo_fixer_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_topo_fixer_SRCS
mris_topo_fixer.cpp
)


add_executable(mris_topo_fixer ${mris_topo_fixer_SRCS})
target_link_libraries(mris_topo_fixer ${FS_LIBS})
install(TARGETS mris_topo_fixer DESTINATION bin)	

