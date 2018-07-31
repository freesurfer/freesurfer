project(mri_concat)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_concat mri_concat.c)
target_link_libraries(mri_concat utils)
install(TARGETS mri_concat DESTINATION bin)
