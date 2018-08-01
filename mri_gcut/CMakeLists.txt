project(mri_gcut)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_gcut mri_gcut.cpp)
add_help(mri_gcut mri_gcut.help.xml)
target_link_libraries(mri_gcut utils)

install(TARGETS mri_gcut DESTINATION bin)
