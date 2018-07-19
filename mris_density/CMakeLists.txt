project(mris_density)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_density mris_density.c)
target_link_libraries(mris_density utils)

install(TARGETS mris_density DESTINATION bin)
