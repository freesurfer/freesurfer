project(mris_volmask)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_volmask mris_volmask.cpp cmd_line_interface.cpp)
target_compile_definitions(mris_volmask PUBLIC NO_VTK)
add_help(mris_volmask mris_volmask.help.xml)
target_link_libraries(mris_volmask utils)
install(TARGETS mris_volmask DESTINATION bin)

add_test_script(NAME mris_volmask_test SCRIPT test_mris_volmask DEPENDS mris_volmask)

add_executable(mris_volmask_novtk mris_volmask_old.cpp cmd_line_interface.cpp)
target_compile_definitions(mris_volmask_novtk PUBLIC NO_VTK)
target_link_libraries(mris_volmask_novtk utils)
install(TARGETS mris_volmask_novtk DESTINATION bin)

if(VTK_FOUND)

  include_directories(${VTK_INCLUDE_DIRS})

  add_executable(mris_volmask_vtk mris_volmask_old.cpp cmd_line_interface.cpp)
  target_link_libraries(mris_volmask_vtk utils ${VTK_LIBRARIES})

  set(WRAPCODE "\
#!/bin/tcsh -ef\n\
source $FREESURFER_HOME/bin/vtk_setup\n\
mris_volmask_vtk.bin $argv"
  )
  install_wrapped(TARGETS mris_volmask_vtk DESTINATION bin WRAPCODE ${WRAPCODE})

endif()
