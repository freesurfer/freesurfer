project(mri_surf2surf)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_surf2surf mri_surf2surf.cpp)
target_link_libraries(mri_surf2surf utils)

add_test_script(NAME mri_surf2surf_test SCRIPT test.sh)

install(TARGETS mri_surf2surf DESTINATION bin)
