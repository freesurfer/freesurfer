project(mris_ms_refine)
include_directories(${mris_ms_refine_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_ms_refine_SRCS
mris_ms_refine.c
)


add_executable(mris_ms_refine ${mris_ms_refine_SRCS})
target_link_libraries(mris_ms_refine ${FS_LIBS})
install(TARGETS mris_ms_refine DESTINATION bin)	

