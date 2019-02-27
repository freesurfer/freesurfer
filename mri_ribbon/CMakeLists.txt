project(mri_ribbon)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_ribbon mri_ribbon.cpp)
target_link_libraries(mri_ribbon utils)

install(TARGETS mri_ribbon DESTINATION bin)
