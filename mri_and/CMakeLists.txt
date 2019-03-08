project(mri_and)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_and mri_and.cpp)
target_link_libraries(mri_and utils)
install(TARGETS mri_and DESTINATION bin)

add_executable(mri_or mri_or.cpp)
target_link_libraries(mri_or utils)
install(TARGETS mri_or DESTINATION bin)
