project(mri_strip_nonwhite)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_strip_nonwhite mri_strip_nonwhite.c)
target_link_libraries(mri_strip_nonwhite utils)

install(TARGETS mri_strip_nonwhite DESTINATION bin)
