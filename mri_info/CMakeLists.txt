project(mri_info)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_info mri_info.cpp)
target_link_libraries(mri_info utils)

install(TARGETS mri_info DESTINATION bin)
