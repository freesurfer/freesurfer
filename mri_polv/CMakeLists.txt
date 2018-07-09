project(mri_polv)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_polv mri_polv.c)
target_link_libraries(mri_polv utils)

install(TARGETS mri_polv DESTINATION bin)
