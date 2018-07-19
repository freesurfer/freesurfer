project(mris_thickness_comparison)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_thickness_comparison mris_thickness_comparison.c)
target_link_libraries(mris_thickness_comparison utils)

install(TARGETS mris_thickness_comparison DESTINATION bin)
