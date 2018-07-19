project(mri_average)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_average mri_average.c)
target_link_libraries(mri_average utils)

install(TARGETS mri_average DESTINATION bin)
