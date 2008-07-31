project(mris_thickness)
include_directories(${mris_thickness_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_thickness_SRCS
mris_cluster_profiles.c
mris_intensity_profile.c
mris_thickness.c
)


add_executable(mris_thickness ${mris_thickness_SRCS})
target_link_libraries(mris_thickness ${FS_LIBS})
install(TARGETS mris_thickness DESTINATION bin)	

