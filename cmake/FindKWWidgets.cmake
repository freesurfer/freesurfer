# KWWidgets Find Module

if(NOT KWWidgets_DIR)
  set(KWWidgets_DIR ${FS_PACKAGES_DIR}/KWWidgets/CVS)
endif()

find_path(KWWidgets_INCLUDE_DIR HINTS ${KWWidgets_DIR} NAMES vtkKWApplication.h PATH_SUFFIXES include/KWWidgets)
find_library(KWWidgets_LIBRARIES HINTS ${KWWidgets_DIR} NAMES KWWidgets PATH_SUFFIXES lib/KWWidgets)
find_package_handle_standard_args(KWWidgets DEFAULT_MSG KWWidgets_INCLUDE_DIR KWWidgets_LIBRARIES)

# make sure the shared libs get installed on linux
if(KWWidgets_FOUND AND NOT APPLE)
  file(GLOB KWW_LIBRARIES_TO_INSTALL "${KWWidgets_DIR}/lib/KWWidgets/lib*.so*")
  install(PROGRAMS ${KWW_LIBRARIES_TO_INSTALL} DESTINATION lib)
endif()
