project(mris_transform)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_transform mris_transform.c)
target_link_libraries(mris_transform utils)

install(TARGETS mris_transform DESTINATION bin)
