project(mri_ms_fitparms)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_ms_fitparms mri_ms_fitparms.cpp)
add_help(mri_ms_fitparms mri_ms_fitparms.help.xml)
target_link_libraries(mri_ms_fitparms utils)

install(TARGETS mri_ms_fitparms DESTINATION bin)
