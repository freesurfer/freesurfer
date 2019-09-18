project(mri_linear_register)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_linear_register mri_linear_register.cpp)
target_link_libraries(mri_linear_register utils)

install(TARGETS mri_linear_register DESTINATION bin)
