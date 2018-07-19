project(mri_ms_fitparms)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_ms_fitparms mri_ms_fitparms.c)
add_help(mri_ms_fitparms mri_ms_fitparms.help.xml)
target_link_libraries(mri_ms_fitparms utils)

install(TARGETS mri_ms_fitparms DESTINATION bin)
