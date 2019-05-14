project(mris_simulate_atrophy)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_simulate_atrophy mris_simulate_atrophy.cpp)
target_link_libraries(mris_simulate_atrophy utils)

install(TARGETS mris_simulate_atrophy DESTINATION bin)
