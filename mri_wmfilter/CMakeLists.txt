project(mri_wmfilter)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_wmfilter mri_wmfilter.c)
target_link_libraries(mri_wmfilter utils)

install(TARGETS mri_wmfilter DESTINATION bin)
