project(mri_convert)

# ATH: temporarily turning off this warning
add_compile_options(-Wno-self-assign)

include_directories(${FS_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/packages/dicom)

add_executable(mri_convert mri_convert.cpp)
add_help(mri_convert mri_convert.help.xml)
target_link_libraries(mri_convert utils)

install(TARGETS mri_convert DESTINATION bin)

add_test_script(NAME mri_convert_test SCRIPT test.sh DEPENDS mri_convert)

add_executable(mri_make_uchar mri_make_uchar.cpp)
target_link_libraries(mri_make_uchar utils)
install(TARGETS mri_make_uchar DESTINATION bin)
