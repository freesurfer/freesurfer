project(mris_warp)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_warp mris_warp.cpp)
target_link_libraries(mris_warp utils)

install(TARGETS mris_warp DESTINATION bin)
