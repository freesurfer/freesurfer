project(mris_spherical_average)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_spherical_average mris_spherical_average.cpp)
target_link_libraries(mris_spherical_average utils)

install(TARGETS mris_spherical_average DESTINATION bin)
