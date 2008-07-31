project(mris_longitudinal_surfaces)
include_directories(${mris_longitudinal_surfaces_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_longitudinal_surfaces_SRCS
mris_longitudinal_surfaces.c
)


add_executable(mris_longitudinal_surfaces ${mris_longitudinal_surfaces_SRCS})
target_link_libraries(mris_longitudinal_surfaces ${FS_LIBS})
install(TARGETS mris_longitudinal_surfaces DESTINATION bin)	


