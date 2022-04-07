project(mri_label2label)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_label2label mri_label2label.cpp)
target_link_libraries(mri_label2label utils)

add_test_script(NAME mri_label2label_test SCRIPT test.sh)

install(TARGETS mri_label2label DESTINATION bin)
