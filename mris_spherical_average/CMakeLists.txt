project(mris_spherical_average)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_spherical_average mris_spherical_average.c)
target_link_libraries(mris_spherical_average utils)

install(TARGETS mris_spherical_average DESTINATION bin)
