project(mris_twoclass)
include_directories(${mris_twoclass_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_twoclass_SRCS
mris_twoclass.c
)


add_executable(mris_twoclass ${mris_twoclass_SRCS})
target_link_libraries(mris_twoclass ${FS_LIBS})
install(TARGETS mris_twoclass DESTINATION bin)	


