project(mris_make_surfaces)
include_directories(${mris_make_surfaces_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_make_surfaces_SRCS
mris_exvivo_surfaces.c
mris_make_surfaces.c
mris_mef_surfaces.c
mris_refine_surfaces.c
)


add_executable(mris_make_surfaces ${mris_make_surfaces_SRCS})
target_link_libraries(mris_make_surfaces ${FS_LIBS})
install(TARGETS mris_make_surfaces DESTINATION bin)	

