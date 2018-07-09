project(mris_anatomical_stats)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_anatomical_stats mris_anatomical_stats.c)
add_help(mris_anatomical_stats mris_anatomical_stats.help.xml)
target_link_libraries(mris_anatomical_stats utils)

install(TARGETS mris_anatomical_stats DESTINATION bin)
