project(AntsN4BiasFieldCorrectionFs)

include_directories(${FS_INCLUDE_DIRS})

add_compile_options(-Wno-self-assign-field)

host_os()
if(HOST_OS MATCHES "Ubuntu18" OR HOST_OS MATCHES "CentOS8")
   message(STATUS "AntsN4BiasFieldCorrectionFs test will be skipped on ${HOST_OS}")
else()
   add_test_script(NAME AntsN4BiasFieldCorrectionFs_test SCRIPT test.sh)
endif()

add_executable(AntsN4BiasFieldCorrectionFs AntsN4BiasFieldCorrectionFs.cpp)
add_help(AntsN4BiasFieldCorrectionFs AntsN4BiasFieldCorrectionFs.help.xml)
target_link_libraries(AntsN4BiasFieldCorrectionFs utils)

install(TARGETS AntsN4BiasFieldCorrectionFs DESTINATION bin)
