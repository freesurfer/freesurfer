project(mri_normalize)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_normalize mri_normalize.cpp)
add_help(mri_normalize mri_normalize.help.xml)
target_link_libraries(mri_normalize utils)
install(TARGETS mri_normalize DESTINATION bin)

host_os()
if(HOST_OS MATCHES "CentOS8")
   message(STATUS "mri_normalize test will be skipped on ${HOST_OS}")
else()
   add_test_script(NAME mri_normalize_test SCRIPT test.sh DEPENDS mri_normalize)
endif()

add_executable(mri_long_normalize mri_long_normalize.cpp)
add_help(mri_long_normalize mri_long_normalize.help.xml)
target_link_libraries(mri_long_normalize utils)
install(TARGETS mri_long_normalize DESTINATION bin)
