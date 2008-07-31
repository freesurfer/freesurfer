project(mris_thickness_comparison)
include_directories(${mris_thickness_comparison_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_thickness_comparison_SRCS
mris_thickness_comparison.c
)


add_executable(mris_thickness_comparison ${mris_thickness_comparison_SRCS})
target_link_libraries(mris_thickness_comparison ${FS_LIBS})
install(TARGETS mris_thickness_comparison DESTINATION bin)	

