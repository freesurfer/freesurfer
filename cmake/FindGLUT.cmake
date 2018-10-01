# GLUT Find Module

if(NOT GLUT_DIR)
  set(GLUT_DIR ${FS_PACKAGES_DIR}/glut/3.7)
endif()

find_path(GLUT_INCLUDE_DIR HINTS ${GLUT_DIR} NAMES GL/glut.h PATH_SUFFIXES include)
find_library(GLUT_LIBRARIES HINTS ${GLUT_DIR} NAMES libglut.a PATH_SUFFIXES lib)
find_package_handle_standard_args(GLUT DEFAULT_MSG GLUT_INCLUDE_DIR GLUT_LIBRARIES)
