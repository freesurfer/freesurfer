project(mri_ribbon)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_ribbon mri_ribbon.c)
target_link_libraries(mri_ribbon utils)

install(TARGETS mri_ribbon DESTINATION bin)
