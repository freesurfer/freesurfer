project(mris_register_label_map)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_register_label_map mris_register_label_map.cpp)
target_link_libraries(mris_register_label_map utils)

install(TARGETS mris_register_label_map DESTINATION bin)
