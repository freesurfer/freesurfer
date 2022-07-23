project(mri_ca_normalize)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_ca_normalize mri_ca_normalize.cpp)
add_help(mri_ca_normalize mri_ca_normalize.help.xml)
target_link_libraries(mri_ca_normalize utils)
install(TARGETS mri_ca_normalize DESTINATION bin)

add_test_script(NAME mri_ca_normalize_test SCRIPT test.sh)
