project(mri_remove_neck)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_remove_neck mri_remove_neck.cpp)
add_help(mri_remove_neck mri_remove_neck.help.xml)
target_link_libraries(mri_remove_neck utils)

install(TARGETS mri_remove_neck DESTINATION bin)
