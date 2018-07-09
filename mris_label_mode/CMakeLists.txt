project(mris_label_mode)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_label_mode mris_label_mode.c)
target_link_libraries(mris_label_mode utils)

install(TARGETS mris_label_mode DESTINATION bin)
