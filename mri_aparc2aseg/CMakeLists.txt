project(mri_aparc2aseg)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_aparc2aseg mri_aparc2aseg.c)
add_help(mri_aparc2aseg mri_aparc2aseg.help.xml)
target_link_libraries(mri_aparc2aseg utils)

install(TARGETS mri_aparc2aseg DESTINATION bin)
