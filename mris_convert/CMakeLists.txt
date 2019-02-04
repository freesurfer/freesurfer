project(mris_convert)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_convert mris_convert.c)
target_link_libraries(mris_convert utils)

add_test_script(NAME mris_convert_test SCRIPT test.sh DEPENDS mris_convert)

install(TARGETS mris_convert DESTINATION bin)
