# Qt Find Module

# default search path
if(USE_QT4)
  set(PACKAGES_QT ${FS_PACKAGES_DIR}/qt/4.8.5)
  if(NOT QT_QMAKE_EXECUTABLE AND EXISTS ${PACKAGES_QT})
    set(QT_INSTALL_DIR ${PACKAGES_QT})
    set(QT_QMAKE_EXECUTABLE ${QT_INSTALL_DIR}/bin/qmake)
  endif()
else()
  set(PACKAGES_QT ${FS_PACKAGES_DIR}/qt/5.6/lib/cmake/Qt5)
  if(NOT Qt5_DIR AND EXISTS ${PACKAGES_QT})
    set(Qt5_DIR ${PACKAGES_QT})
  endif()
  get_filename_component(QT_INSTALL_DIR "${Qt5_DIR}/../../.." ABSOLUTE)
endif()

# find Qt components and trace back root of install
if(USE_QT4)
  set(QT_COMPONENTS QtCore QtGui)
  find_package(Qt4 COMPONENTS ${QT_COMPONENTS})
else()
  set(QT_COMPONENTS Core Widgets)
  if(NOT APPLE)
    set(QT_COMPONENTS ${QT_COMPONENTS} X11Extras)
  endif()
  find_package(Qt5 COMPONENTS ${QT_COMPONENTS})
endif()

# install the shared libraries to the freesurfer lib directory
if(((NOT USE_QT4 AND Qt5_FOUND) OR (USE_QT4 AND Qt4_FOUND)) AND NOT APPLE)
  file(GLOB QT_LIBRARIES_TO_INSTALL "${QT_INSTALL_DIR}/lib/lib*.so*")
  if(QT_LIBRARIES_TO_INSTALL)
    install(PROGRAMS ${QT_LIBRARIES_TO_INSTALL} DESTINATION lib/qt/lib)
    # add Qt library directory to rpath
    set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib/qt/lib:${CMAKE_INSTALL_RPATH}")
    # install the platform plugins as well, and make sure executables know
    #where to find the plugins directory
    if(NOT USE_QT4 AND EXISTS ${QT_INSTALL_DIR}/plugins/platforms)
      install(DIRECTORY ${QT_INSTALL_DIR}/plugins/platforms DESTINATION lib/qt/plugins)
      install(FILES ${CMAKE_SOURCE_DIR}/cmake/qt.conf DESTINATION bin)
    endif()
  endif()
endif()

