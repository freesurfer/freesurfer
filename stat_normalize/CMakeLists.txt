project(stat_normalize)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(stat_normalize stat_normalize.c)
target_link_libraries(stat_normalize utils)

install(TARGETS stat_normalize DESTINATION bin)
