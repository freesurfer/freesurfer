project(mris_label_area)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_label_area mris_label_area.cpp)
target_link_libraries(mris_label_area utils)

install(TARGETS mris_label_area DESTINATION bin)
