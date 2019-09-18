project(mris_map_cuts)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_map_cuts mris_map_cuts.cpp)
target_link_libraries(mris_map_cuts utils)

install(TARGETS mris_map_cuts DESTINATION bin)
