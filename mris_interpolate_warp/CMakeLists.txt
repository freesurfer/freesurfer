project(mris_interpolate_warp)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_interpolate_warp mris_interpolate_warp.c)
target_link_libraries(mris_interpolate_warp utils)

install(TARGETS mris_interpolate_warp DESTINATION bin)
