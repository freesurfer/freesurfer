project(mri_watershed)

include_directories(${FS_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/packages/dicom)

add_executable(mri_watershed mri_watershed.cpp)
add_help(mri_watershed mri_watershed.help.xml)
target_link_libraries(mri_watershed utils)
install(TARGETS mri_watershed DESTINATION bin)

add_test_script(NAME mri_watershed_test SCRIPT test.sh DEPENDS mri_watershed)

add_executable(mri_validate_skull_stripped mri_validate_skull_stripped.cpp)
target_link_libraries(mri_validate_skull_stripped utils)
install(TARGETS mri_validate_skull_stripped DESTINATION bin)

add_subdirectory(brain_volume)
