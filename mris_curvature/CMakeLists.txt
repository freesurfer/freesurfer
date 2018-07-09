project(mris_curvature)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_curvature mris_curvature.c)
add_help(mris_curvature mris_curvature.help.xml)
target_link_libraries(mris_curvature utils)

install(TARGETS mris_curvature DESTINATION bin)
