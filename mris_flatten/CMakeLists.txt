project(mris_flatten)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_flatten mris_flatten.c)
target_link_libraries(mris_flatten utils)
install(TARGETS mris_flatten DESTINATION bin)

