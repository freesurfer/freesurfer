project(mri_mc)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_mc mri_mc.c)
target_link_libraries(mri_mc utils)
install(TARGETS mri_mc DESTINATION bin)

add_executable(mri_pretess mri_pretess.c)
add_help(mri_pretess mri_pretess.help.xml)
target_link_libraries(mri_pretess utils)
install(TARGETS mri_pretess DESTINATION bin)

add_executable(mris_extract_main_component mris_extract_main_component.c)
add_help(mris_extract_main_component mris_extract_main_component.help.xml)
target_link_libraries(mris_extract_main_component utils)
install(TARGETS mris_extract_main_component DESTINATION bin)