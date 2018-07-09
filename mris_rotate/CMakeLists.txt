project(mris_rotate)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_rotate mris_rotate.c)
target_link_libraries(mris_rotate utils)

install(TARGETS mris_rotate DESTINATION bin)
