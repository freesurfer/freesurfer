project(mris_ms_surface_CNR)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_ms_surface_CNR mris_ms_surface_CNR.c)
target_link_libraries(mris_ms_surface_CNR utils)

install(TARGETS mris_ms_surface_CNR DESTINATION bin)
