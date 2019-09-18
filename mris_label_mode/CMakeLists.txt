project(mris_label_mode)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_label_mode mris_label_mode.cpp)
target_link_libraries(mris_label_mode utils)

install(TARGETS mris_label_mode DESTINATION bin)
