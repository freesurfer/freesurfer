project(mri_binarize)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_binarize mri_binarize.cpp)
target_link_libraries(mri_binarize utils)

add_test_script(NAME mri_binarize_test SCRIPT test.sh)

install(TARGETS mri_binarize DESTINATION bin)
