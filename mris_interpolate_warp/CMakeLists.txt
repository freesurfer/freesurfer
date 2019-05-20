project(mris_interpolate_warp)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_interpolate_warp mris_interpolate_warp.cpp)
target_link_libraries(mris_interpolate_warp utils)

install(TARGETS mris_interpolate_warp DESTINATION bin)
