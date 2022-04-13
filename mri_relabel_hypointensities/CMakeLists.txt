project(mri_relabel_hypointensities)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_relabel_hypointensities mri_relabel_hypointensities.cpp)
target_link_libraries(mri_relabel_hypointensities utils)

add_test_script(NAME mri_relabel_hypointensities_test SCRIPT test.sh)

install(TARGETS mri_relabel_hypointensities DESTINATION bin)
