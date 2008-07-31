project(mris_errors)
include_directories(${mris_errors_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_errors_SRCS
mris_errors.c
)


add_executable(mris_errors ${mris_errors_SRCS})
target_link_libraries(mris_errors ${FS_LIBS})
install(TARGETS mris_errors DESTINATION bin)	

