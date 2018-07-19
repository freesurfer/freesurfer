project(mris_simulate_atrophy)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_simulate_atrophy mris_simulate_atrophy.c)
target_link_libraries(mris_simulate_atrophy utils)

install(TARGETS mris_simulate_atrophy DESTINATION bin)
