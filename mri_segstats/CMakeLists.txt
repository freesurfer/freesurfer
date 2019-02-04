project(mri_segstats)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_segstats mri_segstats.c)
add_help(mri_segstats mri_segstats.help.xml)
target_link_libraries(mri_segstats utils)

add_test_script(NAME mri_segstats_test SCRIPT test.sh DEPENDS mri_segstats)

install(TARGETS mri_segstats DESTINATION bin)
