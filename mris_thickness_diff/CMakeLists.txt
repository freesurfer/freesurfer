project(mris_thickness_diff)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_thickness_diff mris_thickness_diff.c)
target_link_libraries(mris_thickness_diff utils)

install(TARGETS mris_thickness_diff DESTINATION bin)
