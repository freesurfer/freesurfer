project(mris_compute_parc_overlap)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_compute_parc_overlap mris_compute_parc_overlap.cpp)
target_link_libraries(mris_compute_parc_overlap utils)

install(TARGETS mris_compute_parc_overlap DESTINATION bin)
