project(mri_aparc2aseg)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_aparc2aseg mri_aparc2aseg.cpp)
add_help(mri_aparc2aseg mri_aparc2aseg.help.xml)
target_link_libraries(mri_aparc2aseg utils)

add_test_script(NAME mri_aparc2aseg_test SCRIPT test.sh)

install(TARGETS mri_aparc2aseg DESTINATION bin)


add_executable(mri_surf2volseg mri_surf2volseg.cpp)
add_help(mri_surf2volseg mri_surf2volseg.help.xml)
target_link_libraries(mri_surf2volseg utils)

#add_test_script(NAME mri_surf2volseg_test SCRIPT test.sh)

install(TARGETS mri_surf2volseg DESTINATION bin)
