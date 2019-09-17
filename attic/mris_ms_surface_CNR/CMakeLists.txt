project(mris_ms_surface_CNR)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_ms_surface_CNR mris_ms_surface_CNR.cpp)
target_link_libraries(mris_ms_surface_CNR utils)

install(TARGETS mris_ms_surface_CNR DESTINATION bin)
