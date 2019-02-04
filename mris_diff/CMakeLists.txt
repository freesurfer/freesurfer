project(mris_diff)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_diff mris_diff.c)
target_link_libraries(mris_diff utils)

add_test_script(NAME mris_diff_test SCRIPT test.sh DEPENDS mris_diff)

install(TARGETS mris_diff DESTINATION bin)
