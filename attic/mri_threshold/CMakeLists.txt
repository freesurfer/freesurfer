project(mri_threshold)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_threshold mri_threshold.cpp)
target_link_libraries(mri_threshold utils)

install(TARGETS mri_threshold DESTINATION bin)
