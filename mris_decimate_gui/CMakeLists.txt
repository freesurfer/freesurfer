project(mris_decimate_gui)

if(BUILD_GUIS AND GTS_FOUND AND VTK_FOUND AND WXWidgets_FOUND)

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${WXWidgets_FLAGS} -Wno-deprecated -Wno-literal-suffix")

  include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/res
    ${FS_INCLUDE_DIRS}
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

  add_executable(mris_decimate_gui.bin ${SOURCES})
  target_link_libraries(mris_decimate_gui.bin
    utils
    fsgdf
    ${VTK_LIBRARIES}
    ${WXWidgets_LIBRARIES}
    ${GTS_LIBRARIES}
    ${OPENGL_LIBRARIES}
    ${X11_LIBRARIES}
  )

  install(PROGRAMS mris_decimate_gui_wrapper DESTINATION bin RENAME mris_decimate_gui)
  install(TARGETS mris_decimate_gui.bin DESTINATION bin) 

endif()
