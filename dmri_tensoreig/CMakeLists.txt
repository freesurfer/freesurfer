project(dmri_tensoreig)
include_directories(${dmri_tensoreig_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(dmri_tensoreig_SRCS
dmri_tensoreig.c
)


add_executable(dmri_tensoreig ${dmri_tensoreig_SRCS})
target_link_libraries(dmri_tensoreig ${FS_LIBS})
install(TARGETS dmri_tensoreig DESTINATION bin)	
