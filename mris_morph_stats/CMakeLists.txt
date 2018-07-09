project(mris_morph_stats)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_morph_stats mris_morph_stats.c)
target_link_libraries(mris_morph_stats utils)

install(TARGETS mris_morph_stats DESTINATION bin)
