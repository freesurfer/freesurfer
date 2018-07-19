project(mri_seghead)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_seghead mri_seghead.c)
target_link_libraries(mri_seghead utils)

install(TARGETS mri_seghead DESTINATION bin)
