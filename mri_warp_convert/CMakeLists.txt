project(mri_warp_convert)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_warp_convert mri_warp_convert.cpp)
add_help(mri_warp_convert mri_warp_convert.help.xml)
target_link_libraries(mri_warp_convert utils)

add_test_script(NAME mri_warp_convert_test SCRIPT test.sh)

install(TARGETS mri_warp_convert DESTINATION bin)
