project(fsgdf)

if(TCLTKTIXBLT_FOUND)

  include_directories(
    ${FS_INCLUDE_DIRS}
    ${MINC_INCLUDE_DIR}
    ${TCLTKTIXBLT_INCLUDE_DIR}
    ${X11_INCLUDE_DIR}
  )

  SET(SOURCES
    fsgdf.c
    fsgdf_wrap.c
    FsgdfPlot.cxx
  )

  add_library(fsgdf STATIC ${SOURCES})

  target_link_libraries(fsgdf utils ${MINC_LIBRARIES} ${NetCDF_LIBRARIES})

endif()
