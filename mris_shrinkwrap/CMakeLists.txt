project(mris_shrinkwrap)
include_directories(${mris_shrinkwrap_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_shrinkwrap_SRCS
mris_AA_shrinkwrap.c
mris_shrinkwrap.c
)


add_executable(mris_shrinkwrap ${mris_shrinkwrap_SRCS})
target_link_libraries(mris_shrinkwrap ${FS_LIBS})
install(TARGETS mris_shrinkwrap DESTINATION bin)	


