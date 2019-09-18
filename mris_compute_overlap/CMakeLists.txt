project(mris_compute_overlap)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_compute_overlap mris_compute_overlap.cpp)
target_link_libraries(mris_compute_overlap utils)

install(TARGETS mris_compute_overlap DESTINATION bin)
