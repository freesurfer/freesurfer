project(mris_smooth)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_smooth mris_smooth.cpp)
add_help(mris_smooth mris_smooth.help.xml)
target_link_libraries(mris_smooth utils)

add_test_script(NAME mris_smooth_posix_test SCRIPT test.sh)

install(TARGETS mris_smooth DESTINATION bin)
