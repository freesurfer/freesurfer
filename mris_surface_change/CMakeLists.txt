project(mris_surface_change)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_surface_change mris_surface_change.cpp)
target_link_libraries(mris_surface_change utils)

install(TARGETS mris_surface_change DESTINATION bin)
