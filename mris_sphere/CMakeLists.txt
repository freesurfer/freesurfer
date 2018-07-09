project(mris_sphere)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_sphere mris_sphere.c)
add_help(mris_sphere mris_sphere.help.xml)
target_link_libraries(mris_sphere utils)
install(TARGETS mris_sphere DESTINATION bin)

add_executable(mris_remove_negative_vertices mris_remove_negative_vertices.c)
target_link_libraries(mris_remove_negative_vertices utils)
install(TARGETS mris_remove_negative_vertices DESTINATION bin)
