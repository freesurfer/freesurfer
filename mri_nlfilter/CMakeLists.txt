project(mri_nlfilter)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_nlfilter mri_nlfilter.c)
target_link_libraries(mri_nlfilter utils)

install(TARGETS mri_nlfilter DESTINATION bin)
