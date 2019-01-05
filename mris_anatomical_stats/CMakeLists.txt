project(mris_anatomical_stats)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_anatomical_stats mris_anatomical_stats.c)
add_help(mris_anatomical_stats mris_anatomical_stats.help.xml)
target_link_libraries(mris_anatomical_stats utils)

add_test_script(NAME mris_anatomical_stats_test SCRIPT test_mris_anatomical_stats DEPENDS mris_anatomical_stats)
configure_file(testdata2.tar.gz . COPYONLY)
configure_file(test1.expected . COPYONLY)
configure_file(test2.expected . COPYONLY)
configure_file(test2_logfile.expected . COPYONLY)
configure_file(test3.expected . COPYONLY)
configure_file(test6.expected . COPYONLY)

install(TARGETS mris_anatomical_stats DESTINATION bin)
