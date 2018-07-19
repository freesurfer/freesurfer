project(mri_cc)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_cc mri_cc.c)
add_help(mri_cc mri_cc.help.xml)
target_link_libraries(mri_cc utils)

install(TARGETS mri_cc DESTINATION bin)
