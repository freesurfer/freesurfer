project(mris_decimate_gui)

if(BUILD_GUIS AND GTS_FOUND AND VTK_FOUND AND WXWidgets_FOUND)

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${WXWidgets_FLAGS} -Wno-deprecated -Wno-literal-suffix")

  include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/res
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/mris_decimate
    ${WXWidgets_INCLUDE_DIRS}
    ${GTS_INCLUDE_DIRS}
    ${VTK_INCLUDE_DIRS})

  set(SOURCES
    ../mris_decimate/mris_decimate.cpp
    DecimatePanel.cpp
    MainApp.cpp
    MainWindow.cpp
    MyCmdLineParser.cpp
    RenderPanel.cpp
    wxVTKRenderWindowInteractor.cpp
    res/MainWindowBase.cpp
    res/DecimatePanelBase.cpp
  )

  add_executable(mris_decimate_gui ${SOURCES})
  target_link_libraries(mris_decimate_gui
    utils
    fsgdf
    ${VTK_LIBRARIES}
    ${WXWidgets_LIBRARIES}
    ${GTS_LIBRARIES}
    ${OPENGL_LIBRARIES} X11 Xmu
  )

  set(WRAPCODE "\
#!/bin/tcsh -ef\n\
source $FREESURFER_HOME/bin/vtk_setup\n\
if (-e /bin/pwd) then\n\
  setenv HOME `cd $HOME && /bin/pwd`\n\
  setenv SUBJECTS_DIR `cd $SUBJECTS_DIR && /bin/pwd`\n\
endif\n\
mris_decimate_gui.bin $argv"
  )
  install_wrapped(TARGETS mris_decimate_gui DESTINATION bin WRAPCODE ${WRAPCODE})

endif()
