project(mris_show)
include_directories(${mris_show_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_show_SRCS
mris_show.c
)


add_executable(mris_show ${mris_show_SRCS})
target_link_libraries(mris_show ${FS_LIBS})
install(TARGETS mris_show DESTINATION bin)	

