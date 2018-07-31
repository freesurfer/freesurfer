# KWWidgets Find Module

if(NOT KWWidgets_DIR)
  set(KWWidgets_DIR ${FS_PACKAGES_DIR}/KWWidgets/CVS)
endif()

find_path(KWWidgets_INCLUDE_DIR HINTS ${KWWidgets_DIR} NAMES vtkKWApplication.h PATH_SUFFIXES include/KWWidgets)
find_library(KWWidgets_LIBRARIES HINTS ${KWWidgets_DIR} NAMES KWWidgets PATH_SUFFIXES lib/KWWidgets)
find_package_handle_standard_args(KWWidgets DEFAULT_MSG KWWidgets_INCLUDE_DIR KWWidgets_LIBRARIES)

# make sure the kwwidgets path gets linked from the lib directory during install
symlink(${KWWidgets_DIR} ${CMAKE_INSTALL_PREFIX}/lib/KWWidgets)
