project(mri_compute_overlap)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_compute_overlap mri_compute_overlap.c)
target_link_libraries(mri_compute_overlap utils)

install(TARGETS mri_compute_overlap DESTINATION bin)
