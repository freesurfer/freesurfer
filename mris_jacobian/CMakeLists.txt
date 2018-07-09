project(mris_jacobian)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_jacobian mris_jacobian.c)
add_help(mris_jacobian mris_jacobian.help.xml)
target_link_libraries(mris_jacobian utils)

install(TARGETS mris_jacobian DESTINATION bin)
