project(mri_compile_edits)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_compile_edits mri_compile_edits.cpp)
target_link_libraries(mri_compile_edits utils)

install(TARGETS mri_compile_edits DESTINATION bin)
