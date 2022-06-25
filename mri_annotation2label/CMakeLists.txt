project(mri_annotation2label)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_annotation2label mri_annotation2label.cpp)
target_link_libraries(mri_annotation2label utils)

add_test_script(NAME mri_annotation2label_test SCRIPT test.sh)

install(TARGETS mri_annotation2label DESTINATION bin)
