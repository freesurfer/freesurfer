project(mris_thickness_comparison)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_thickness_comparison mris_thickness_comparison.cpp)
target_link_libraries(mris_thickness_comparison utils)

install(TARGETS mris_thickness_comparison DESTINATION bin)
