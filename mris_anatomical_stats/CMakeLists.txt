project(mris_anatomical_stats)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_anatomical_stats mris_anatomical_stats.cpp)
add_help(mris_anatomical_stats mris_anatomical_stats.help.xml)
target_link_libraries(mris_anatomical_stats utils)

add_test_script(NAME mris_anatomical_stats_test SCRIPT test.sh)

install(TARGETS mris_anatomical_stats DESTINATION bin)
