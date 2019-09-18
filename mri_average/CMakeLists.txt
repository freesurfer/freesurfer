project(mri_average)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_average mri_average.cpp)
target_link_libraries(mri_average utils)

install(TARGETS mri_average DESTINATION bin)
