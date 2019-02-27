project(mri_simulate_atrophy)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_simulate_atrophy mri_simulate_atrophy.cpp)
target_link_libraries(mri_simulate_atrophy utils)

install(TARGETS mri_simulate_atrophy DESTINATION bin)
