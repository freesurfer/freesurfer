project(mris_convert)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_convert mris_convert.c)
target_link_libraries(mris_convert utils)

add_test_script(NAME mris_convert_test SCRIPT test_mris_convert DEPENDS mris_convert)
# todo: test gifti

install(TARGETS mris_convert DESTINATION bin)
