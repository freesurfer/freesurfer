project(mri_threshold)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_threshold mri_threshold.c)
target_link_libraries(mri_threshold utils)

install(TARGETS mri_threshold DESTINATION bin)
