project(mris_distance_to_label)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_distance_to_label mris_distance_to_label.cpp)
target_link_libraries(mris_distance_to_label utilscpp)
install(TARGETS mris_distance_to_label DESTINATION bin)
