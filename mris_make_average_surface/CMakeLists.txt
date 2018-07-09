project(mris_make_average_surface)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_make_average_surface mris_make_average_surface.c)
target_link_libraries(mris_make_average_surface utils)

install(TARGETS mris_make_average_surface DESTINATION bin)
