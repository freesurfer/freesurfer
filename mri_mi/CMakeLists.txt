project(mri_mi)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_mi mri_mi.c)
target_link_libraries(mri_mi utils)

install(TARGETS mri_mi DESTINATION bin)
