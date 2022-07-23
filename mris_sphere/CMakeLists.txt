project(mris_sphere)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_sphere mris_sphere.cpp)
add_help(mris_sphere mris_sphere.help.xml)
target_link_libraries(mris_sphere utils)
install(TARGETS mris_sphere DESTINATION bin)

add_test_script(NAME mris_sphere_test SCRIPT test.sh)

add_executable(mris_remove_negative_vertices mris_remove_negative_vertices.cpp)
target_link_libraries(mris_remove_negative_vertices utils)
install(TARGETS mris_remove_negative_vertices DESTINATION bin)
