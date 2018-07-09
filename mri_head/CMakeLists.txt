project(mri_head)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_head mri_head.c)
target_link_libraries(mri_head utils)

install(TARGETS mri_head DESTINATION bin)
