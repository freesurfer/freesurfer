project(mris_entropy)
include_directories(${mris_entropy_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_entropy_SRCS
mris_entropy.c
)


add_executable(mris_entropy ${mris_entropy_SRCS})
target_link_libraries(mris_entropy ${FS_LIBS})
install(TARGETS mris_entropy DESTINATION bin)	

