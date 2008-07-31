project(.-test_makevol)
include_directories(${.-test_makevol_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(.-test_makevol_SRCS
main.c
)


add_executable(.-test_makevol ${.-test_makevol_SRCS})
target_link_libraries(.-test_makevol ${FS_LIBS})
install(TARGETS .-test_makevol DESTINATION bin)	

