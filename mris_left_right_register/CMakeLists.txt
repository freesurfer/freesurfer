project(mris_left_right_register)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_left_right_register mris_left_right_register.c)
target_link_libraries(mris_left_right_register utils)

install(TARGETS mris_left_right_register DESTINATION bin)
