project(mri_cc)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_cc mri_cc.cpp)
add_help(mri_cc mri_cc.help.xml)
target_link_libraries(mri_cc utils)

install(TARGETS mri_cc DESTINATION bin)
