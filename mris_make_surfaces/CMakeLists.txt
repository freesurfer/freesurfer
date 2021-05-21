project(mris_make_surfaces)

add_compile_options(-Wno-deprecated -Wno-overloaded-virtual -Wno-inconsistent-missing-override -Wno-unused-variable -Wno-unused-local-typedefs -Wno-reorder)

include_directories(${FS_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/resurf/Code)

add_executable(mris_make_surfaces mris_make_surfaces.cpp)
add_help(mris_make_surfaces mris_make_surfaces.help.xml)
target_link_libraries(mris_make_surfaces utils)
install(TARGETS mris_make_surfaces DESTINATION bin)

host_os()
if(NOT APPLE)
   if(HOST_OS MATCHES "CentOS8")
      message(STATUS "mris_make_surfaces test will be skipped on ${HOST_OS}")
      message(STATUS "mris_place_surface test will be skipped on ${HOST_OS}")
   else()
      add_test_script(NAME mris_make_surfaces_test SCRIPT test.sh DEPENDS mris_make_surfaces)
      add_test_script(NAME mris_place_surface_test SCRIPT test.sh DEPENDS mris_place_surface)
      # commands in the same subdir that can share the same input/output files should not run in ||
      set_tests_properties(mris_make_surfaces_test mris_place_surface_test PROPERTIES RUN_SERIAL TRUE)
   endif()
# test commands exit with non-zero status from too many vertex differences on MacOS
else()
   message(STATUS "mris_make_surfaces test will be skipped on MacOS")
   message(STATUS "mris_place_surface test will be skipped on MacOS")
endif()

add_executable(mris_place_surface mris_place_surface.cpp)
add_help(mris_place_surface mris_place_surface.help.xml)
target_link_libraries(mris_place_surface utils)
install(TARGETS mris_place_surface DESTINATION bin)

add_executable(mris_autodet_gwstats mris_autodet_gwstats.cpp)
target_link_libraries(mris_autodet_gwstats utils)
install(TARGETS mris_autodet_gwstats DESTINATION bin)

add_executable(mris_refine_surfaces mris_refine_surfaces.cpp)
target_link_libraries(mris_refine_surfaces utils)
install(TARGETS mris_refine_surfaces DESTINATION bin)

add_executable(mris_mef_surfaces mris_mef_surfaces.cpp)
target_link_libraries(mris_mef_surfaces utils)
install(TARGETS mris_mef_surfaces DESTINATION bin)

if(NOT MINIMAL)
   add_executable(mris_exvivo_surfaces mris_exvivo_surfaces.cpp)
   target_link_libraries(mris_exvivo_surfaces utils)
   install(TARGETS mris_exvivo_surfaces DESTINATION bin)
endif()
