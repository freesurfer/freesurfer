project(mri_warp_convert)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_warp_convert mri_warp_convert.cpp)
add_help(mri_warp_convert mri_warp_convert.help.xml)
target_link_libraries(mri_warp_convert utils)

# 2/2024 - fails on arm64 MacOS builds
if(NOT APPLE_ARM64)
   add_test_script(NAME mri_warp_convert_test SCRIPT test.sh)
endif()

install(TARGETS mri_warp_convert DESTINATION bin)
