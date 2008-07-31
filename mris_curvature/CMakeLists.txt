project(mris_curvature)
include_directories(${mris_curvature_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_curvature_SRCS
mris_curvature.c
)


add_executable(mris_curvature ${mris_curvature_SRCS})
target_link_libraries(mris_curvature ${FS_LIBS})
install(TARGETS mris_curvature DESTINATION bin)	

