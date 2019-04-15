project(mri_compute_overlap)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_compute_overlap mri_compute_overlap.cpp)
target_link_libraries(mri_compute_overlap utils)

install(TARGETS mri_compute_overlap DESTINATION bin)
