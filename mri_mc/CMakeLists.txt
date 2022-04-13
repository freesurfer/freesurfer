project(mri_mc)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_mc mri_mc.cpp)
target_link_libraries(mri_mc utils)
install(TARGETS mri_mc DESTINATION bin)

add_test_script(NAME mri_mc_test SCRIPT test.sh)

add_executable(mri_pretess mri_pretess.cpp)
add_help(mri_pretess mri_pretess.help.xml)
target_link_libraries(mri_pretess utils)
install(TARGETS mri_pretess DESTINATION bin)

add_executable(mris_extract_main_component mris_extract_main_component.cpp)
add_help(mris_extract_main_component mris_extract_main_component.help.xml)
target_link_libraries(mris_extract_main_component utils)
install(TARGETS mris_extract_main_component DESTINATION bin)
