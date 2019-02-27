project(mris_sample_label)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_sample_label mris_sample_label.cpp)
target_link_libraries(mris_sample_label utils)

install(TARGETS mris_sample_label DESTINATION bin)
