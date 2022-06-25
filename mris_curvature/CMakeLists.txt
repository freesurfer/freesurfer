project(mris_curvature)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_curvature mris_curvature.cpp)
add_help(mris_curvature mris_curvature.help.xml)
target_link_libraries(mris_curvature utils)

add_test_script(NAME mris_curvature_test SCRIPT test.sh)

install(TARGETS mris_curvature DESTINATION bin)
