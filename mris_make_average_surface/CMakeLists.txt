project(mris_make_average_surface)
include_directories(${mris_make_average_surface_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_make_average_surface_SRCS
mris_make_average_surface.c
)


add_executable(mris_make_average_surface ${mris_make_average_surface_SRCS})
target_link_libraries(mris_make_average_surface ${FS_LIBS})
install(TARGETS mris_make_average_surface DESTINATION bin)	

