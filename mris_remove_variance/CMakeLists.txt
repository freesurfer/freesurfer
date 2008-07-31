project(mris_remove_variance)
include_directories(${mris_remove_variance_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_remove_variance_SRCS
mris_remove_variance.c
)


add_executable(mris_remove_variance ${mris_remove_variance_SRCS})
target_link_libraries(mris_remove_variance ${FS_LIBS})
install(TARGETS mris_remove_variance DESTINATION bin)	

