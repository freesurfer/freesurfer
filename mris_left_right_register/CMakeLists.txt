project(mris_left_right_register)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_left_right_register mris_left_right_register.cpp)
target_link_libraries(mris_left_right_register utils)
install(TARGETS mris_left_right_register DESTINATION bin)
