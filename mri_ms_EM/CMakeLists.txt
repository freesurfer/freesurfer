project(mri_ms_EM)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_ms_EM mri_ms_EM.c PoissonSolver.c)
target_link_libraries(mri_ms_EM utils)
install(TARGETS mri_ms_EM DESTINATION bin)

