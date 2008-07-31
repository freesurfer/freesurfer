project(mris_classify_thickness)
include_directories(${mris_classify_thickness_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_classify_thickness_SRCS
mris_classify_thickness.c
)


add_executable(mris_classify_thickness ${mris_classify_thickness_SRCS})
target_link_libraries(mris_classify_thickness ${FS_LIBS})
install(TARGETS mris_classify_thickness DESTINATION bin)	

