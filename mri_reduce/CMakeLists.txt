project(mri_reduce)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_reduce mri_reduce.cpp)
target_link_libraries(mri_reduce utils)

install(TARGETS mri_reduce DESTINATION bin)
