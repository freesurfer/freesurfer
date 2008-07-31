project(.-test_window_env)
include_directories(${.-test_window_env_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(.-test_window_env_SRCS
test_window_env.c
)


add_executable(.-test_window_env ${.-test_window_env_SRCS})
target_link_libraries(.-test_window_env ${FS_LIBS})
install(TARGETS .-test_window_env DESTINATION bin)	

