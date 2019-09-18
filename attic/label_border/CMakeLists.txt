project(label_border)

include_directories(${FS_INCLUDE_DIRS})

add_executable(label_border label_border.cpp)
target_link_libraries(label_border utils)

install(TARGETS label_border DESTINATION bin)
