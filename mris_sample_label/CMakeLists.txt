project(mris_sample_label)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_sample_label mris_sample_label.c)
target_link_libraries(mris_sample_label utils)

install(TARGETS mris_sample_label DESTINATION bin)
