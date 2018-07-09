project(mri_histo_eq)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_histo_eq mri_histo_eq.c)
target_link_libraries(mri_histo_eq utils)

install(TARGETS mri_histo_eq DESTINATION bin)
