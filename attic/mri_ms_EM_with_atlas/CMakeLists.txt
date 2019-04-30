project(mri_ms_EM_with_atlas)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_ms_EM_with_atlas mri_ms_EM_with_atlas.cpp PoissonSolver.cpp)
target_link_libraries(mri_ms_EM_with_atlas utils)
install(TARGETS mri_ms_EM_with_atlas DESTINATION bin)
