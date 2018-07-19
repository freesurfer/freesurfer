project(mri_update_gca)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_update_gca mri_update_gca.c)
target_link_libraries(mri_update_gca utils)

install(TARGETS mri_update_gca DESTINATION bin)
