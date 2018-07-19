project(mri_reduce)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_reduce mri_reduce.c)
target_link_libraries(mri_reduce utils)

install(TARGETS mri_reduce DESTINATION bin)
