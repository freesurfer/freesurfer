project(mris_make_average_surface)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_make_average_surface mris_make_average_surface.cpp)
target_link_libraries(mris_make_average_surface utils)

install(TARGETS mris_make_average_surface DESTINATION bin)
