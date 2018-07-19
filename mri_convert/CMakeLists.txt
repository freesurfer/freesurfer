project(mri_convert)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_convert mri_convert.c)
add_help(mri_convert mri_convert.help.xml)
target_link_libraries(mri_convert utils fsgdf)
install_wrapped(TARGETS mri_convert DESTINATION bin)

add_test_script(NAME mri_convert_test SCRIPT test_mri_convert DEPENDS mri_convert)

add_executable(mri_make_uchar mri_make_uchar.c)
target_link_libraries(mri_make_uchar utils)
install(TARGETS mri_make_uchar DESTINATION bin)
