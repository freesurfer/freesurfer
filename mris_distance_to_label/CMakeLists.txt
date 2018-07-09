project(mris_distance_to_label)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_distance_to_label mris_distance_to_label.cpp)
target_link_libraries(mris_distance_to_label utils)

install(TARGETS mris_distance_to_label DESTINATION bin)
