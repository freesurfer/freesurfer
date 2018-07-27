project(mri_convert)

include_directories(${FS_INCLUDE_DIRS} ${MINC_INCLUDE_DIR})

add_executable(mri_convert mri_convert.c)
add_help(mri_convert mri_convert.help.xml)
target_link_libraries(mri_convert utils fsgdf)

install(TARGETS mri_convert DESTINATION bin)

add_test_script(NAME mri_convert_test SCRIPT test_mri_convert DEPENDS mri_convert)

add_executable(mri_make_uchar mri_make_uchar.c)
target_link_libraries(mri_make_uchar utils)
install(TARGETS mri_make_uchar DESTINATION bin)
