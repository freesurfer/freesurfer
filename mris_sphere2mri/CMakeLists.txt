project(mris_make_surfaces)

add_compile_options(-Wno-deprecated -Wno-overloaded-virtual -Wno-inconsistent-missing-override -Wno-unused-variable -Wno-unused-local-typedefs -Wno-reorder)

include_directories(${FS_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/resurf/Code)

add_executable(mris_sphere2mri mris_sphere2mri.cpp)
add_help(mris_sphere2mri mris_sphere2mri.help.xml mri_stopmask.help.xml)
target_link_libraries(mris_sphere2mri utils)
install(TARGETS mris_sphere2mri DESTINATION bin)

host_os()
if(NOT APPLE)
   if(HOST_OS MATCHES "CentOS8")
      message(STATUS "Skipping mris_sphere2mri test on ${HOST_OS}")
   else()
      add_test_script(NAME mris_sphere2mri_test SCRIPT test.sh)
      # commands in the same subdir that can share the same input/output files should not run in ||
      set_tests_properties(mris_sphere2mri_test PROPERTIES RUN_SERIAL TRUE)
   endif()
# test commands exit with non-zero status from too many vertex differences on MacOS
else()
   message(STATUS "Skipping mris_sphere2mri test on MacOS")
endif()

