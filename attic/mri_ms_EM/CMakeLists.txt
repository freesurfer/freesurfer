project(mri_ms_EM)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_ms_EM mri_ms_EM.cpp PoissonSolver.cpp)
target_link_libraries(mri_ms_EM utils)
install(TARGETS mri_ms_EM DESTINATION bin)

