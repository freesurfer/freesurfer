project(mris_fill)
include_directories(${mris_fill_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_fill_SRCS
mris_fill.c
)


add_executable(mris_fill ${mris_fill_SRCS})
target_link_libraries(mris_fill ${FS_LIBS})
install(TARGETS mris_fill DESTINATION bin)	

