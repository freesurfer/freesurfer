project(mris_curvature)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_curvature mris_curvature.c)
add_help(mris_curvature mris_curvature.help.xml)
target_link_libraries(mris_curvature utils)

add_test_script(NAME mris_curvature_test SCRIPT test_mris_curvature DEPENDS mris_curvature)

install(TARGETS mris_curvature DESTINATION bin)
