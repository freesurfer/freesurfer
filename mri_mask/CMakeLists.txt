project(mri_mask)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_mask mri_mask.cpp)
add_help(mri_mask mri_mask.help.xml)
target_link_libraries(mri_mask utils)

add_test_script(NAME mri_mask_test SCRIPT test.sh DEPENDS mri_mask)
install(TARGETS mri_mask DESTINATION bin)

