project(mris_decimate_gui)

if(GTS_FOUND AND VTK_FOUND AND WXWidgets_FOUND)

  add_definitions(-D__WXGTK__ -D__WXGTK20__)

  include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/res
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/mris_decimate
    ${WXWIDGETS_INCLUDE_DIRS}
    ${GTS_INCLUDE_DIRS}
    ${VTK_INCLUDE_DIR})

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
  target_link_libraries(mris_decimate_gui utils fsgdf ${VTK_LIBRARIES} ${WXWidgets_LIBRARIES} ${GTS_LIBRARIES} ${OPENGL_LIBRARIES} X11 Xmu SM Xxf86vm
  gtk-x11-2.0
  gdk-x11-2.0
  atk-1.0
  gdk_pixbuf-2.0
  pangoxft-1.0
  pangox-1.0
  pango-1.0
  gobject-2.0
  gmodule-2.0
  gthread-2.0
  Xinerama
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
