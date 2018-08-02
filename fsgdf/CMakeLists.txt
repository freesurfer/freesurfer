project(fsgdf)

include_directories(
  ${FS_INCLUDE_DIRS}
  ${MINC_INCLUDE_DIR}
  ${X11_INCLUDE_DIR}
)

if(TCLTKTIXBLT_FOUND)
  include_directories(${TCLTKTIXBLT_INCLUDE_DIR})
  set(SOURCES fsgdf_wrap.c FsgdfPlot.cxx)
endif()

add_library(fsgdf STATIC fsgdf.c ${SOURCES})
target_link_libraries(fsgdf utils ${MINC_LIBRARIES} ${NetCDF_LIBRARIES})

if(TCLTKTIXBLT_FOUND AND NOT APPLE)
  add_executable(test_fsgdf EXCLUDE_FROM_ALL test.c fsgdf_wrap.c)
  target_link_libraries(test_fsgdf fsgdf ${TCLTKTIXBLT_LIBRARIES})
  add_test_script(NAME fsgdf_test SCRIPT test.py DEPENDS test_fsgdf)
endif()

