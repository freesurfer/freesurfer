project(mris_fill)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_fill mris_fill.c)
target_link_libraries(mris_fill utils)

install(TARGETS mris_fill DESTINATION bin)
