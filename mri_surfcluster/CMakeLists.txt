project(mri_surfcluster)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_surfcluster mri_surfcluster.c)
target_link_libraries(mri_surfcluster utils)
install_wrapped(TARGETS mri_surfcluster DESTINATION bin)
