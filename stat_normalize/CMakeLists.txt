project(stat_normalize)

include_directories(${FS_INCLUDE_DIRS})

add_executable(stat_normalize stat_normalize.cpp)
target_link_libraries(stat_normalize utils)

install(TARGETS stat_normalize DESTINATION bin)
