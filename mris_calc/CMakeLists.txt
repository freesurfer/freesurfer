project(mris_calc)
include_directories(${mris_calc_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_calc_SRCS
mris_calc.c
)


add_executable(mris_calc ${mris_calc_SRCS})
target_link_libraries(mris_calc ${FS_LIBS})
install(TARGETS mris_calc DESTINATION bin)	
