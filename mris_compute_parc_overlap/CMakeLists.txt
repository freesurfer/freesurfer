project(mris_compute_parc_overlap)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_compute_parc_overlap mris_compute_parc_overlap.c)
target_link_libraries(mris_compute_parc_overlap utils)

install(TARGETS mris_compute_parc_overlap DESTINATION bin)
