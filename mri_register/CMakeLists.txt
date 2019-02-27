project(mri_register)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_register mri_register.cpp)
target_link_libraries(mri_register utils)

install(TARGETS mri_register DESTINATION bin)
