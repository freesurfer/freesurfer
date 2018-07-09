project(stim_polar)

include_directories(${CMAKE_SOURCE_DIR}/include ${GLUT_INCLUDE_DIR}/GL)

add_executable(stim_polar stim_polar.c)
target_link_libraries(stim_polar utils)

install(TARGETS stim_polar DESTINATION bin)
