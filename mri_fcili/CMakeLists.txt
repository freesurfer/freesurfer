project(mri_fcili)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_fcili mri_fcili.c)
target_link_libraries(mri_fcili utils)

install(TARGETS mri_fcili DESTINATION bin)
