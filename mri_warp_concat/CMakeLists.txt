project(mri_warp_concat)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_warp_concat mri_warp_concat.cpp)
add_help(mri_warp_concat mri_warp_concat.help.xml)
target_link_libraries(mri_warp_concat utils)

install(TARGETS mri_warp_concat DESTINATION bin)
