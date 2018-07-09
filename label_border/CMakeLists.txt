project(label_border)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(label_border label_border.c)
target_link_libraries(label_border utils)

install(TARGETS label_border DESTINATION bin)
