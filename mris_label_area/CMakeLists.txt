project(mris_label_area)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_label_area mris_label_area.c)
target_link_libraries(mris_label_area utils)

install(TARGETS mris_label_area DESTINATION bin)
