project(mri_compile_edits)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_compile_edits mri_compile_edits.c)
target_link_libraries(mri_compile_edits utils)

install(TARGETS mri_compile_edits DESTINATION bin)
