project(mris_transform)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_transform mris_transform.c)
target_link_libraries(mris_transform utils)

add_test_script(NAME mris_transform_test SCRIPT test DEPENDS mris_transform)

install(TARGETS mris_transform DESTINATION bin)
