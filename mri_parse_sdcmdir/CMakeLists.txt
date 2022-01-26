project(mri_parse_sdcmdir)

include_directories(${FS_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/packages/dicom ${CMAKE_SOURCE_DIR}/packages/dcm2niix)

add_executable(mri_parse_sdcmdir mri_parse_sdcmdir.cpp)
target_link_libraries(mri_parse_sdcmdir utils)

add_test_script(NAME mri_parse_sdcmdir_test SCRIPT test.sh DEPENDS mri_parse_sdcmdir)

install(TARGETS mri_parse_sdcmdir DESTINATION bin)
