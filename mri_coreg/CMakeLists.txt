project(mri_coreg)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_coreg mri_coreg.cpp)
target_link_libraries(mri_coreg utils)

add_test_script(NAME mri_coreg_test SCRIPT test.sh)

install(TARGETS mri_coreg DESTINATION bin)
