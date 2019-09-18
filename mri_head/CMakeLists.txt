project(mri_head)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_head mri_head.cpp)
target_link_libraries(mri_head utils)

install(TARGETS mri_head DESTINATION bin)
