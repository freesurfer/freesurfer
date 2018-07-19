project(mri_paint)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_paint mri_paint.c)
target_link_libraries(mri_paint utils)

install(TARGETS mri_paint DESTINATION bin)
