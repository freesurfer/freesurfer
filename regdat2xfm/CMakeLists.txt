project(regdat2xfm)
include_directories(${regdat2xfm_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(regdat2xfm_SRCS
regdat2xfm.cpp
)


add_executable(regdat2xfm ${regdat2xfm_SRCS})
target_link_libraries(regdat2xfm ${FS_LIBS})
install(TARGETS regdat2xfm DESTINATION bin)	

