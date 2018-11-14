# Qt Find Module

# default search path
set(PACKAGES_QT ${FS_PACKAGES_DIR}/qt/5.6/lib/cmake/Qt5)
if(NOT Qt5_DIR AND EXISTS ${PACKAGES_QT})
  set(Qt5_DIR ${PACKAGES_QT})
endif()

# set Qt components to search for
set(QT_COMPONENTS Core Widgets)
if(NOT APPLE)
  set(QT_COMPONENTS ${QT_COMPONENTS} X11Extras)
endif()

find_package(Qt5 COMPONENTS ${QT_COMPONENTS})
if(Qt5_FOUND)
  # trace back the root of the qt installation
  get_filename_component(QT_INSTALL_DIR "${Qt5_DIR}/../../.." ABSOLUTE)

  if(NOT APPLE)
    # install all the shared libraries to the freesurfer lib directory
    file(GLOB QT_LIBRARIES_TO_INSTALL "${QT_INSTALL_DIR}/lib/lib*.so*")
    install(PROGRAMS ${QT_LIBRARIES_TO_INSTALL} DESTINATION lib/qt/lib)
    # add Qt library directory to rpath
    set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib/qt/lib:${CMAKE_INSTALL_RPATH}")
    # install the platform plugins as well
    install(DIRECTORY ${QT_INSTALL_DIR}/plugins/platforms DESTINATION lib/qt/plugins)
    # make sure executables know where to find our plugin directory
    file(WRITE ${CMAKE_INSTALL_PREFIX}/bin/qt.conf "[Paths]\nPlugins=../lib/qt/plugins\n")
  endif()

endif()
