project(unix)

include_directories(
  ${CMAKE_SOURCE_DIR}/include
  ${MINC_INCLUDE_DIR}
)

SET(SOURCES
  cuserid.c
  ftime.c
  glStub.c
  oglutil.c
  texfont.c
  unix_thread.c
  windiag.c
  xvmri.c
  xvutil.c
  xwin.c
  xwindow.c
)

add_library(unix STATIC ${SOURCES})

target_link_libraries(unix ${MINC_LIBRARIES})
add_external_dependencies(unix minc)