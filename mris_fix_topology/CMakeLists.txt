project(mris_fix_topology)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_fix_topology mris_fix_topology.cpp)
add_help(mris_fix_topology mris_fix_topology.help.xml)
target_link_libraries(mris_fix_topology utils)
install(TARGETS mris_fix_topology DESTINATION bin)

host_os()
if(HOST_OS MATCHES "CentOS8")
   message(STATUS "mris_fix_topology test will be skipped on ${HOST_OS}")
else()
   add_test_script(NAME mris_fix_topology_test SCRIPT test.sh DEPENDS mris_fix_topology)
endif()

add_executable(mris_remove_intersection mris_remove_intersection.cpp)
add_help(mris_remove_intersection mris_remove_intersection.help.xml)
target_link_libraries(mris_remove_intersection utils)
install(TARGETS mris_remove_intersection DESTINATION bin)
