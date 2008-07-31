project(mris_average_curvature)
include_directories(${mris_average_curvature_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_average_curvature_SRCS
mris_average_curvature.c
)


add_executable(mris_average_curvature ${mris_average_curvature_SRCS})
target_link_libraries(mris_average_curvature ${FS_LIBS})
install(TARGETS mris_average_curvature DESTINATION bin)	


