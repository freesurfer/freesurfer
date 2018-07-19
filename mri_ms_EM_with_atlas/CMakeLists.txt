project(mri_ms_EM_with_atlas)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_ms_EM_with_atlas mri_ms_EM_with_atlas.c PoissonSolver.c)
target_link_libraries(mri_ms_EM_with_atlas utils)
install(TARGETS mri_ms_EM_with_atlas DESTINATION bin)
