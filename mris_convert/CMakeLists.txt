project(mris_convert)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_convert mris_convert.c)
target_link_libraries(mris_convert utils)

add_test_script(NAME mris_convert_test SCRIPT test_mris_convert DEPENDS mris_convert)
add_test_script(NAME mris_convert_gifti_test SCRIPT test_gifti DEPENDS mris_convert)

install(TARGETS mris_convert DESTINATION bin)
