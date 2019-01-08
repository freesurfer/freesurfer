# Qt Find Module

if(NOT Qt5_DIR)
  # default search path
  if(EXISTS ${FS_PACKAGES_DIR}/qt/5.11/lib/cmake/Qt5)
    set(Qt5_DIR ${FS_PACKAGES_DIR}/qt/5.11/lib/cmake/Qt5)
  elseif(EXISTS ${FS_PACKAGES_DIR}/qt/5.6/lib/cmake/Qt5)
    set(Qt5_DIR ${FS_PACKAGES_DIR}/qt/5.6/lib/cmake/Qt5)
  endif()
endif()
get_filename_component(QT_INSTALL_DIR "${Qt5_DIR}/../../.." ABSOLUTE)

# find Qt components and trace back root of install
set(QT_COMPONENTS Core Widgets)
if(NOT APPLE)
  set(QT_COMPONENTS ${QT_COMPONENTS} X11Extras)
endif()

find_package(Qt5 COMPONENTS ${QT_COMPONENTS})

# install the shared libraries to the freesurfer lib directory
if(Qt5_FOUND AND NOT APPLE)
  file(GLOB QT_LIBRARIES_TO_INSTALL "${QT_INSTALL_DIR}/lib/lib*.so*")
  if(QT_LIBRARIES_TO_INSTALL)
    install(PROGRAMS ${QT_LIBRARIES_TO_INSTALL} DESTINATION lib/qt/lib)
    # add Qt library directory to rpath
    set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib/qt/lib:${CMAKE_INSTALL_RPATH}")
    # install the platform plugins as well, and make sure executables know
    # where to find the plugins directory
    if(EXISTS ${QT_INSTALL_DIR}/plugins/platforms)
      install(DIRECTORY ${QT_INSTALL_DIR}/plugins/platforms DESTINATION lib/qt/plugins)
      install(FILES ${CMAKE_SOURCE_DIR}/qt/qt.conf DESTINATION bin)
    endif()
  endif()
endif()
