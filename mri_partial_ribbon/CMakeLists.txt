project(mri_partial_ribbon)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_partial_ribbon mri_partial_ribbon.cpp)
target_link_libraries(mri_partial_ribbon utils)

install(TARGETS mri_partial_ribbon DESTINATION bin)
