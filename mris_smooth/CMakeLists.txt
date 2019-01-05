project(mris_smooth)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_smooth mris_smooth.c)
add_help(mris_smooth mris_smooth.help.xml)
target_link_libraries(mris_smooth utils)

add_test_script(NAME mris_smooth_posix_test SCRIPT test_fs_posix DEPENDS mris_smooth)
configure_file(testdata.tar.gz . COPYONLY)

install(TARGETS mris_smooth DESTINATION bin)
