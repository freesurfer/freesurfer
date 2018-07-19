project(mris_euler_number)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_euler_number mris_euler_number.c)
add_help(mris_euler_number mris_euler_number.help.xml)
target_link_libraries(mris_euler_number utils)

install(TARGETS mris_euler_number DESTINATION bin)
