project(mri_surf2vol)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_surf2vol mri_surf2vol.cpp)
target_link_libraries(mri_surf2vol utils)

add_test_script(NAME mri_surf2vol_test SCRIPT test.sh)

install(TARGETS mri_surf2vol DESTINATION bin)
