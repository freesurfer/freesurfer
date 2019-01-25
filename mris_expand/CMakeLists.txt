project(mris_expand)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_expand mris_expand.c)
target_link_libraries(mris_expand utils)

add_test_script(NAME mris_expand_test SCRIPT test.sh DEPENDS mris_expand)

install(TARGETS mris_expand DESTINATION bin)
