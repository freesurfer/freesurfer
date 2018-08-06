project(unix)

include_directories(
  ${FS_INCLUDE_DIRS}
  ${MINC_INCLUDE_DIR}
  ${X11_INCLUDE_DIR}
)

set(SOURCES
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

target_link_libraries(unix ${MINC_LIBRARIES} ${NetCDF_LIBRARIES} ${OPENGL_LIBRARIES})
