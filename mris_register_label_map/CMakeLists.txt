project(mris_register_label_map)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_register_label_map mris_register_label_map.c)
target_link_libraries(mris_register_label_map utils)

install(TARGETS mris_register_label_map DESTINATION bin)
