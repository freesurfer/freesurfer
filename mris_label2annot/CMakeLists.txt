project(mris_label2annot)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_label2annot mris_label2annot.cpp)
target_link_libraries(mris_label2annot utils)

add_test_script(NAME mris_label2annot_test SCRIPT test.sh)

install(TARGETS mris_label2annot DESTINATION bin)
