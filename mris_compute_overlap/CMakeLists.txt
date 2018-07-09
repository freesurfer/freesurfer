project(mris_compute_overlap)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_compute_overlap mris_compute_overlap.c)
target_link_libraries(mris_compute_overlap utils)

install(TARGETS mris_compute_overlap DESTINATION bin)
