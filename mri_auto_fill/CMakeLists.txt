project(mri_auto_fill)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_auto_fill mri_auto_fill.c)
target_link_libraries(mri_auto_fill utils)

install(TARGETS mri_auto_fill DESTINATION bin)
