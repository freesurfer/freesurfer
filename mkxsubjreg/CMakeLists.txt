project(mkxsubjreg)
include_directories(${mkxsubjreg_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mkxsubjreg_SRCS
mkxsubjreg.c
)


add_executable(mkxsubjreg ${mkxsubjreg_SRCS})
target_link_libraries(mkxsubjreg ${FS_LIBS})
install(TARGETS mkxsubjreg DESTINATION bin)	

