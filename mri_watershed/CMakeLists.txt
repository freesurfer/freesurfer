project(mri_watershed)

include_directories(${FS_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/packages/dicom ${CMAKE_SOURCE_DIR}/packages/dcm2niix)

add_executable(mri_watershed mri_watershed.cpp)
add_help(mri_watershed mri_watershed.help.xml)
target_link_libraries(mri_watershed utils)
install(TARGETS mri_watershed DESTINATION bin)

host_os()
if(HOST_OS MATCHES "CentOS8")
   message(STATUS "mri_watershed test will be skipped on ${HOST_OS}")
else()
   add_test_script(NAME mri_watershed_test SCRIPT test.sh)
endif()

add_executable(mri_validate_skull_stripped mri_validate_skull_stripped.cpp)
target_link_libraries(mri_validate_skull_stripped utils)
install(TARGETS mri_validate_skull_stripped DESTINATION bin)

add_subdirectory(brain_volume)
