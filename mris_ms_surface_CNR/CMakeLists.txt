project(mris_ms_surface_CNR)
include_directories(${mris_ms_surface_CNR_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_ms_surface_CNR_SRCS
mris_ms_surface_CNR.c
)


add_executable(mris_ms_surface_CNR ${mris_ms_surface_CNR_SRCS})
target_link_libraries(mris_ms_surface_CNR ${FS_LIBS})
install(TARGETS mris_ms_surface_CNR DESTINATION bin)	


