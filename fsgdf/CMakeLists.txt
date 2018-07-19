project(fsgdf)

include_directories(
  ${CMAKE_SOURCE_DIR}/include
  ${MINC_INCLUDE_DIR}
)

SET(SOURCES
  fsgdf.c
  fsgdf_wrap.c
  FsgdfPlot.cxx
)

add_library(fsgdf STATIC ${SOURCES})

target_link_libraries(fsgdf utils ${MINC_LIBRARIES} ${NetCDF_LIBRARIES})
