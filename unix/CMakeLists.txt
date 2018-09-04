project(unix)

include_directories(
  ${FS_INCLUDE_DIRS}
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

target_link_libraries(unix minc netcdf ${OPENGL_LIBRARIES})
