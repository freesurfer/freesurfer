project(mris_diff)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_diff mris_diff.c)
target_link_libraries(mris_diff utils)

install(TARGETS mris_diff DESTINATION bin)
