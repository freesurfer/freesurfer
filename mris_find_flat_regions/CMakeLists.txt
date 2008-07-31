project(mris_find_flat_regions)
include_directories(${mris_find_flat_regions_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_find_flat_regions_SRCS
mris_find_flat_regions.c
)


add_executable(mris_find_flat_regions ${mris_find_flat_regions_SRCS})
target_link_libraries(mris_find_flat_regions ${FS_LIBS})
install(TARGETS mris_find_flat_regions DESTINATION bin)	

