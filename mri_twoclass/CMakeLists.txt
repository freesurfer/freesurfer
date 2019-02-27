project(mri_twoclass)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_twoclass mri_twoclass.cpp)
target_link_libraries(mri_twoclass utils)

install(TARGETS mri_twoclass DESTINATION bin)
