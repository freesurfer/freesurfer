project(tksurfer)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(tksurfer tksurfer.c)
target_link_libraries(tksurfer utils)

install(TARGETS tksurfer DESTINATION bin)
