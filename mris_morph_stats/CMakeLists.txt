project(mris_morph_stats)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_morph_stats mris_morph_stats.cpp)
target_link_libraries(mris_morph_stats utils)

install(TARGETS mris_morph_stats DESTINATION bin)
