project(mris_warp)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_warp mris_warp.c)
target_link_libraries(mris_warp utils)

install(TARGETS mris_warp DESTINATION bin)
