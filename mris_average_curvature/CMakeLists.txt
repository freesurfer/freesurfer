project(mris_average_curvature)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_average_curvature mris_average_curvature.c)
target_link_libraries(mris_average_curvature utils)

install(TARGETS mris_average_curvature DESTINATION bin)
