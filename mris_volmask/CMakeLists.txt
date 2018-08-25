project(mris_volmask)

include_directories(${FS_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")

add_executable(mris_volmask mris_volmask.cpp cmd_line_interface.cpp)
target_compile_definitions(mris_volmask PRIVATE NO_VTK)
add_help(mris_volmask mris_volmask.help.xml)
target_link_libraries(mris_volmask utils)
install(TARGETS mris_volmask DESTINATION bin)

add_test_script(NAME mris_volmask_test SCRIPT test_mris_volmask DEPENDS mris_volmask)

add_executable(mris_volmask_novtk mris_volmask_old.cpp cmd_line_interface.cpp)
target_compile_definitions(mris_volmask_novtk PRIVATE NO_VTK)
add_help(mris_volmask_novtk mris_volmask.help.xml)
target_link_libraries(mris_volmask_novtk utils)
install(TARGETS mris_volmask_novtk DESTINATION bin)

if(VTK_FOUND AND NOT MINIMAL)

  include_directories(${VTK_INCLUDE_DIRS})

  add_executable(mris_volmask_vtk.bin mris_volmask_old.cpp cmd_line_interface.cpp)
  add_help(mris_volmask_vtk.bin mris_volmask.help.xml)
  target_link_libraries(mris_volmask_vtk.bin utils ${VTK_LIBRARIES})

  install(PROGRAMS mris_volmask_vtk_wrapper DESTINATION bin RENAME mris_volmask_vtk)
  install(TARGETS mris_volmask_vtk.bin DESTINATION bin)

endif()
