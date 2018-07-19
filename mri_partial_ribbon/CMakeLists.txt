project(mri_partial_ribbon)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_partial_ribbon mri_partial_ribbon.c)
target_link_libraries(mri_partial_ribbon utils)

install(TARGETS mri_partial_ribbon DESTINATION bin)
